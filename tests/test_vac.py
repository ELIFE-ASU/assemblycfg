import json
from collections import Counter
from pathlib import Path

import pytest

import assemblycfg as cfg
from assemblycfg import vac


@pytest.fixture
def no_vac(monkeypatch):
    def unavailable(*args, **kwargs):
        raise cfg.VacError("vac disabled for this test")

    monkeypatch.setattr(vac, "solve_vac", unavailable)


@pytest.fixture(scope="session")
def real_vac(tmp_path_factory):
    """A working vac, installed into a temporary cache if not already present."""
    cache = tmp_path_factory.mktemp("cache")
    patch = pytest.MonkeyPatch()
    patch.setenv("XDG_CACHE_HOME", str(cache))
    try:
        path = cfg.find_vac()
    except cfg.VacError as error:
        patch.undo()
        pytest.skip(f"vac is unavailable: {error}")
    patch.setenv("VAC_PATH", path)
    yield path
    patch.undo()


def closed_form(data, **kwargs):
    return cfg.vac_lower_bound(data, use_vac=False, **kwargs)


def test_mol_unit_counts_include_explicit_hydrogens():
    mol = cfg.smi_to_mol("[H]C#C[H]")
    assert cfg.mol_unit_counts(mol) == Counter({("C", "H", 1): 2, ("C", "C", 3): 1})
    assert cfg.mol_unit_counts(mol, strip_hydrogen=True) == Counter({("C", "C", 3): 1})


def test_mol_unit_counts_leave_the_input_unchanged():
    graph = cfg.smi_to_nx("CCO")
    edges = graph.number_of_edges()
    cfg.mol_unit_counts(graph, strip_hydrogen=True)
    assert graph.number_of_edges() == edges


def test_string_unit_counts():
    assert cfg.string_unit_counts("abracadabra") == Counter(
        {"a": 5, "b": 2, "r": 2, "c": 1, "d": 1})


@pytest.mark.parametrize(
    ("n", "length"),
    [(1, 0), (2, 1), (7, 4), (1000, 12), (9999, 16), (10_000, 14), (2**20, 20)],
)
def test_scalar_chain_length(n, length):
    """Exact through 9999 (A003313), then ceil(log2(n)), which never exceeds l(n)."""
    assert cfg.scalar_chain_length(n) == length


@pytest.mark.parametrize(
    ("data", "expected"),
    [
        pytest.param("a", 0, id="unit"),
        pytest.param("abab", 2, id="abab"),
        pytest.param("a" * 7, 4, id="power"),
        pytest.param("abcd", 3, id="support"),
        pytest.param(["ab", "ba", "aab"], 2, id="joint-distinct-targets"),
    ],
)
def test_closed_form_bound_for_strings(data, expected):
    assert closed_form(data) == expected


def test_disconnected_components_are_separate_targets():
    """Two hydrogen molecules share their only bond, so the joint index is 0."""
    assert closed_form(cfg.smi_to_nx("[H][H].[H][H]")) == 0


def test_merging_rare_units_never_raises_the_bound():
    units = list(range(40))
    vectors = [[i % 5 + 1 for i in units]]
    merged_units, merged = vac._merge_rare_units(units, vectors, 32)
    assert len(merged_units) == len(merged[0]) == 32
    assert sum(merged[0]) == sum(vectors[0])
    assert vac._closed_form_bound(merged) <= vac._closed_form_bound(vectors)


def test_unproven_chain_length_is_never_returned(monkeypatch):
    """An unproven chain is an upper bound on chain length, not on assembly index."""
    record = {"length": 20, "lower_bound": 5, "proven_optimal": False}
    monkeypatch.setattr(vac, "solve_vac", lambda vectors, timeout: record)
    bound, info = cfg.vac_lower_bound("abcabd", return_info=True)
    assert bound == 5
    assert info["source"] == "vac" and not info["proven_optimal"]

    record = {"length": 5, "lower_bound": 4, "proven_optimal": True}
    bound, info = cfg.vac_lower_bound("abcabd", return_info=True)
    assert bound == 5 and info["proven_optimal"]


def test_missing_vac_falls_back_with_a_warning(no_vac):
    with pytest.warns(UserWarning, match="vac is unavailable"):
        bound, info = cfg.vac_lower_bound("abab", return_info=True)
    assert bound == 2 and info["source"] == "closed_form"


def test_entries_beyond_vac_limits_skip_vac(no_vac):
    with pytest.warns(UserWarning, match="exceed vac's limit"):
        assert cfg.vac_lower_bound("ab" * 1100) == 14  # l(2200)


@pytest.mark.parametrize("data", [42, [], ["ab", cfg.smi_to_nx("CC")]])
def test_rejects_unsupported_and_mixed_input(data):
    with pytest.raises(TypeError):
        cfg.vac_lower_bound(data, use_vac=False)


@pytest.fixture
def fake_cargo(tmp_path, monkeypatch):
    """Stand in for cargo, recording each install command."""
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    monkeypatch.delenv("VAC_REF", raising=False)
    monkeypatch.setattr(vac.shutil, "which", lambda name: f"/usr/bin/{name}")
    calls = []

    class Done:
        def __init__(self, returncode):
            self.returncode, self.stderr = returncode, "no such ref"

    def run(command, **kwargs):
        calls.append(command)
        if "--branch" in command and command[command.index("--branch") + 1] == "v1":
            return Done(101)  # v1 is a tag, not a branch
        root = Path(command[command.index("--root") + 1])
        (root / "bin").mkdir(parents=True, exist_ok=True)
        (root / "bin" / vac._VAC_EXECUTABLE).write_text("vac")
        return Done(0)

    monkeypatch.setattr(vac.subprocess, "run", run)
    return calls


def test_install_vac_caches_by_ref(tmp_path, fake_cargo):
    path = cfg.install_vac()
    assert path == str(tmp_path / "assemblycfg" / "vac" / "bin" / vac._VAC_EXECUTABLE)
    assert json.loads((tmp_path / "assemblycfg" / "vac" / "install.json").read_text()) == {
        "repository": vac.VAC_REPOSITORY, "ref": None}
    assert cfg.install_vac() == path
    assert len(fake_cargo) == 1, "a cached install must not reinstall"

    cfg.install_vac(ref="v1")
    assert [c[c.index("--git") + 2:c.index("--root")] for c in fake_cargo[1:]] == [
        ["--branch", "v1"], ["--tag", "v1"]]
    cfg.install_vac(ref="80ba219")
    assert "--rev" in fake_cargo[-1]


def test_install_vac_reports_missing_cargo(tmp_path, monkeypatch):
    monkeypatch.setattr(vac.shutil, "which", lambda name: None)
    monkeypatch.setattr(vac.Path, "home", lambda: tmp_path)
    with pytest.raises(cfg.VacError, match="cargo was not found"):
        cfg.install_vac(force=True)


def test_find_vac_prefers_vac_path(monkeypatch):
    monkeypatch.setenv("VAC_PATH", "/configured/vac")
    assert cfg.find_vac() == "/configured/vac"


@pytest.mark.parametrize("smiles", ["CCO", "c1ccccc1", "OC(=O)CCC(=O)O", "NCC(=O)NCC(=O)O"])
def test_lower_bound_is_below_det_upper_bound(real_vac, smiles):
    graph = cfg.smi_to_nx(smiles)
    upper, _, _ = cfg.calculate_assembly_path_det(graph)
    assert closed_form(graph) <= cfg.vac_lower_bound(graph) <= upper


@pytest.mark.parametrize("data", ["abracadabra", "bbcbaabab", "aaaadbbbbcaa",
                                  ["abcabc", "abcab"]])
def test_lower_bound_is_below_repair_upper_bound(real_vac, data):
    upper, _, _ = cfg.repair_with_pathways(data)
    assert cfg.vac_lower_bound(data) <= upper


def test_vac_proves_tight_bounds(real_vac):
    """Where the chain is the pathway, vac proves the index exactly."""
    bound, info = cfg.vac_lower_bound("abracadabra", return_info=True)
    assert (bound, info["proven_optimal"], info["source"]) == (7, True, "vac")
    bound, info = cfg.vac_lower_bound(cfg.smi_to_nx("CCO"), return_info=True)
    assert (bound, info["proven_optimal"]) == (6, True)
