"""
Lower bounds on assembly index from vector addition chains.

Counting how many copies of each basic unit an object contains maps it to a
non-negative integer vector: bond types (element pair and bond order) for
molecules and characters for strings. The map is a homomorphism from assembly
space to vector addition: a basic unit maps to a unit vector, a join maps to
the sum of its parts' vectors, and reusing an object reuses a chain element.
Every assembly pathway therefore maps to a vector addition chain of the same
length, so the minimal chain length bounds the assembly index from below.
Several objects share one chain, which bounds their joint assembly index.

Chains are solved by ``vac`` from
`additionchains <https://github.com/ELIFE-ASU/additionchains>`_, installed on
demand with ``cargo``. Without it, the closed-form bounds that ``vac`` starts
from are computed here instead.
"""

import json
import os
import platform
import re
import shutil
import subprocess
import warnings
from collections import Counter
from functools import cache
from pathlib import Path
from typing import Any, Dict, Hashable, List, Optional, Sequence, Tuple, Union

import networkx as nx
from rdkit import Chem

from .molecules import mol_to_nx, remove_hydrogen_from_graph

VAC_REPOSITORY = "https://github.com/ELIFE-ASU/additionchains"
_VAC_EXECUTABLE = "vac.exe" if platform.system() == "Windows" else "vac"
# Limits of vac's fixed-width vectors.
_VAC_MAX_DIM = 32
_VAC_MAX_ENTRY = 1024
_CHAIN_TABLE = Path(__file__).parent / "data" / "integer_chain_9999.txt"

Object = Union[str, nx.Graph, Chem.Mol]


class VacError(RuntimeError):
    """``vac`` could not be found, installed or run, or it rejected its input."""


# --- Locating and installing vac -------------------------------------------

def _vac_cache_dir() -> Path:
    """Return ``<cache>/assemblycfg/vac``, honouring ``XDG_CACHE_HOME``."""
    root = os.environ.get("XDG_CACHE_HOME") or "~/.cache"
    return Path(root).expanduser() / "assemblycfg" / "vac"


def _find_cargo() -> str:
    """Find cargo on ``PATH``, then in rustup's default ``~/.cargo/bin``."""
    cargo = shutil.which("cargo")
    if cargo is None:
        suffix = ".exe" if platform.system() == "Windows" else ""
        candidate = Path.home() / ".cargo" / "bin" / ("cargo" + suffix)
        if candidate.is_file() and os.access(candidate, os.X_OK):
            cargo = str(candidate)
    if cargo is None:
        raise VacError("cargo was not found, so vac cannot be installed. Install Rust "
                       "from https://rustup.rs, or set VAC_PATH to a vac executable.")
    return cargo


def _ref_arguments(ref: Optional[str]) -> List[List[str]]:
    """cargo selects a git revision by kind; try each kind *ref* could be."""
    if ref is None:
        return [[]]
    if re.fullmatch(r"[0-9a-f]{7,40}", ref):
        return [["--rev", ref]]
    return [["--branch", ref], ["--tag", ref]]


def install_vac(ref: Optional[str] = None, force: bool = False) -> str:
    """
    Install vac into the assemblycfg cache with ``cargo install``.

    Parameters
    ----------
    ref : str, optional
        Branch, tag or commit of additionchains to install. Defaults to
        ``VAC_REF`` if set, otherwise the repository's default branch.
        Changing the ref reinstalls a cached executable.
    force : bool, optional
        Reinstall even when the cache already holds the requested ref.

    Returns
    -------
    str
        Path to the installed executable,
        ``$XDG_CACHE_HOME/assemblycfg/vac/bin/vac`` (``~/.cache`` by default).

    Raises
    ------
    VacError
        If cargo is missing or the installation fails.
    """
    ref = ref or os.environ.get("VAC_REF") or None
    prefix = _vac_cache_dir()
    executable = prefix / "bin" / _VAC_EXECUTABLE
    record_file = prefix / "install.json"
    record = {"repository": VAC_REPOSITORY, "ref": ref}
    try:
        cached = json.loads(record_file.read_text()) == record
    except (OSError, ValueError):
        cached = False
    if cached and not force and executable.is_file():
        return str(executable)

    cargo = _find_cargo()
    # A rustup cargo finds rustc only when its own directory is on PATH.
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([str(Path(cargo).parent), env.get("PATH", "")])
    errors = []
    for arguments in _ref_arguments(ref):
        command = [cargo, "install", "--locked", "--force", "--quiet",
                   "--git", VAC_REPOSITORY, *arguments,
                   "--root", str(prefix), "vac-cli"]
        print(f"Installing vac from {VAC_REPOSITORY} ({ref or 'default branch'})",
              flush=True)
        done = subprocess.run(command, env=env, capture_output=True, text=True)
        if done.returncode == 0:
            record_file.write_text(json.dumps(record) + "\n")
            return str(executable)
        errors.append(done.stderr.strip())
    raise VacError("Installing vac failed:\n" + "\n".join(errors))


def find_vac() -> str:
    """
    Return the path to the vac executable, installing it if necessary.

    Resolution order is ``VAC_PATH``, ``vac`` on ``PATH``, the cached install,
    and finally :func:`install_vac`.
    """
    return os.environ.get("VAC_PATH") or shutil.which("vac") or install_vac()


def solve_vac(vectors: Sequence[Sequence[int]],
              timeout: float = 10.0,
              vac_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Run ``vac solve`` on target vectors that share one chain.

    Parameters
    ----------
    vectors : Sequence[Sequence[int]]
        Non-negative integer target vectors of equal length.
    timeout : float, optional
        Budget for vac's exact search in seconds; 0 means no limit.
    vac_path : str, optional
        Executable to run. Defaults to :func:`find_vac`.

    Returns
    -------
    dict
        vac's JSON record. ``length`` is the chain found, ``lower_bound`` the
        proven bound it searched up from, and ``proven_optimal`` whether
        ``length`` is minimal.

    Raises
    ------
    VacError
        If vac rejects the input, fails, or overruns its budget.
    """
    command = [vac_path or find_vac(), "solve",
               *(",".join(str(int(c)) for c in v) for v in vectors),
               "--json", "--timeout", str(float(timeout))]
    # vac stops its own search at the timeout; allow a margin before giving up.
    wall = float(timeout) + 10.0 if timeout else None
    try:
        done = subprocess.run(command, capture_output=True, text=True, timeout=wall)
    except subprocess.TimeoutExpired:
        raise VacError(f"vac did not exit within {wall:.0f} seconds") from None
    # 0 is proven minimal and 2 found but unproven; both print a record.
    if done.returncode not in (0, 2):
        raise VacError(done.stderr.strip() or f"vac exited with {done.returncode}")
    return json.loads(done.stdout)


# --- Mapping objects to vectors --------------------------------------------

def _as_graph(mol: Union[nx.Graph, Chem.Mol], strip_hydrogen: bool) -> nx.Graph:
    """Return the coloured molecular graph, leaving the input unchanged."""
    if isinstance(mol, Chem.Mol):
        # Standardization modifies its input; it also adds hydrogens.
        graph = mol_to_nx(Chem.Mol(mol))
    elif isinstance(mol, nx.Graph):
        graph = mol.copy()
    else:
        raise TypeError(f"Input not supported: {type(mol).__name__}")
    return remove_hydrogen_from_graph(graph) if strip_hydrogen else graph


def _graph_unit_counts(graph: nx.Graph) -> Counter:
    """Count edges by (sorted endpoint colours, edge colour)."""
    counts = Counter()
    for u, v, data in graph.edges(data=True):
        a, b = sorted((graph.nodes[u].get("color"), graph.nodes[v].get("color")), key=str)
        counts[(a, b, data.get("color"))] += 1
    return counts


def mol_unit_counts(mol: Union[nx.Graph, Chem.Mol],
                    strip_hydrogen: bool = False) -> Counter:
    """
    Count the bonds of each type in a molecule.

    Parameters
    ----------
    mol : nx.Graph or Chem.Mol
        A molecular graph with node ``color`` atom symbols and edge ``color``
        bond orders, or an RDKit molecule, which is converted with
        :func:`mol_to_nx` and so kekulized and given explicit hydrogens.
    strip_hydrogen : bool, optional
        If True, remove hydrogen atoms first.

    Returns
    -------
    Counter
        Maps ``(symbol, symbol, bond_order)``, symbols sorted, to the number
        of such bonds.

    Examples
    --------
    >>> import assemblycfg as cfg
    >>> cfg.mol_unit_counts(cfg.smi_to_mol("C=C"), strip_hydrogen=True)
    Counter({('C', 'C', 2): 1})
    """
    return _graph_unit_counts(_as_graph(mol, strip_hydrogen))


def string_unit_counts(string: str) -> Counter:
    """Count the occurrences of each character in a string."""
    return Counter(string)


def _targets(data: Union[Object, Sequence[Object]], strip_hydrogen: bool) -> List[Counter]:
    """Convert the input to one unit count per assembly target."""
    if isinstance(data, (str, nx.Graph, Chem.Mol)):
        items = [data]
    elif isinstance(data, Sequence) and data:
        items = list(data)
    else:
        raise TypeError("Input must be a string, graph, molecule or non-empty list of them")

    if all(isinstance(item, str) for item in items):
        return [string_unit_counts(item) for item in items]
    if any(isinstance(item, str) for item in items):
        raise TypeError("Cannot mix strings with graphs or molecules")

    targets = []
    for item in items:
        graph = _as_graph(item, strip_hydrogen)
        # A disconnected molecule's components are separate targets: its
        # joint assembly index need not build their union.
        targets.extend(_graph_unit_counts(graph.subgraph(component))
                       for component in nx.connected_components(graph))
    return targets


def _to_vectors(targets: List[Counter]) -> Tuple[List[Hashable], List[List[int]]]:
    """Lay the targets out as vectors over one shared, deterministic unit order."""
    totals = Counter()
    for target in targets:
        totals.update(target)
    units = sorted(totals, key=lambda unit: (-totals[unit], str(unit)))
    vectors = [[target.get(unit, 0) for unit in units] for target in targets]
    # Zero vectors need no step, and vac rejects them.
    return units, [v for v in vectors if any(v)]


def _merge_rare_units(units: List[Hashable], vectors: List[List[int]],
                      max_dim: int) -> Tuple[List[Hashable], List[List[int]]]:
    """
    Sum the rarest units into one coordinate so at most *max_dim* remain.

    Summing coordinates is a linear map taking unit vectors to unit vectors,
    so it carries any chain to a chain no longer than itself, and the bound
    for the merged vectors still holds for the originals.
    """
    if len(units) <= max_dim:
        return units, vectors
    keep = max_dim - 1  # units are sorted most frequent first
    merged = [v[:keep] + [sum(v[keep:])] for v in vectors]
    return units[:keep] + [tuple(units[keep:])], merged


# --- Closed-form bounds ------------------------------------------------------

@cache
def _chain_table() -> List[int]:
    """Minimal scalar addition chain lengths l(n) for n = 0..9999."""
    table = [0]
    with open(_CHAIN_TABLE) as file:
        for line in file:
            fields = line.split()
            if fields and fields[0].isdigit():
                table.append(int(fields[3]))
    return table


def scalar_chain_length(n: int) -> int:
    """
    Minimal addition chain length l(n), or a lower bound beyond the table.

    Exact for n up to 9999; above that, ``ceil(log2(n))``, which never exceeds
    l(n).
    """
    table = _chain_table()
    return table[n] if n < len(table) else (n - 1).bit_length()


def _closed_form_bound(vectors: List[List[int]]) -> int:
    """
    Combine vac's closed-form lower bounds on a shared vector addition chain.

    Per target ``t`` with nonzero coordinates ``a_1..a_r``:

    - sum: coordinate sums form a scalar chain containing ``sum(t)``.
    - Olivos duality: ``max(max l(a_i), #distinct a_i > 1) + r - 1``, which
      subsumes the projection (``l(a_i)``) and support (``r - 1``) bounds.

    Across targets, each distinct non-unit target needs its own step.
    """
    bound = 0
    distinct = set()
    for vector in vectors:
        coords = [c for c in vector if c]
        total = sum(coords)
        if total > 1:
            distinct.add(tuple(vector))
        sequence = max(max(scalar_chain_length(c) for c in coords),
                       len({c for c in coords if c > 1}))
        bound = max(bound, scalar_chain_length(total), sequence + len(coords) - 1)
    return max(bound, len(distinct))


# --- Public bound ------------------------------------------------------------

def vac_lower_bound(data: Union[Object, Sequence[Object]],
                    strip_hydrogen: bool = False,
                    timeout: float = 10.0,
                    use_vac: bool = True,
                    return_info: bool = False) -> Union[int, Tuple[int, Dict[str, Any]]]:
    """
    Bound the assembly index from below with a vector addition chain.

    The object is mapped to the vector of its unit multiplicities (bond types
    for molecules, characters for strings), and the minimal addition chain
    reaching that vector is a lower bound on its assembly index. A list of
    objects shares one chain, bounding their joint assembly index.

    Parameters
    ----------
    data : str, nx.Graph, Chem.Mol, or a list of them
        One object, or several of the same kind to bound jointly. Each
        connected component of a molecule is a separate target.
    strip_hydrogen : bool, optional
        If True, remove hydrogens from molecules first.
    timeout : float, optional
        Budget in seconds for vac's exact search; 0 means no limit. A search
        that runs out still yields a valid, possibly weaker, bound.
    use_vac : bool, optional
        If False, skip vac and return only the closed-form bound.
    return_info : bool, optional
        If True, also return a dict describing the calculation.

    Returns
    -------
    int or (int, dict)
        The lower bound. With ``return_info``, also a dict with keys ``units``
        (coordinate labels), ``vectors``, ``source`` (``"vac"`` or
        ``"closed_form"``), ``proven_optimal`` (whether the bound is the exact
        minimal chain length) and ``vac`` (vac's record, or None).

    Notes
    -----
    vac's chain length is used only when vac proves it minimal; an unproven
    chain is an upper bound on chain length and says nothing about assembly
    index. Inputs beyond vac's limits are reduced first: units beyond 32 are
    merged, and any entry above 1024 falls back to the closed-form bound, as
    does a missing vac, with a warning.

    Examples
    --------
    >>> import assemblycfg as cfg
    >>> cfg.vac_lower_bound("abab", use_vac=False)
    2
    """
    units, vectors = _to_vectors(_targets(data, strip_hydrogen))
    info = {"units": units, "vectors": vectors, "source": "closed_form",
            "proven_optimal": False, "vac": None}

    bound = _closed_form_bound(vectors) if vectors else 0
    if all(sum(v) <= 1 for v in vectors):
        # Only basic units are targets, and they are free.
        info["proven_optimal"] = True
    elif use_vac:
        merged_units, merged = _merge_rare_units(units, vectors, _VAC_MAX_DIM)
        if max(max(v) for v in merged) > _VAC_MAX_ENTRY:
            warnings.warn(f"Unit counts exceed vac's limit of {_VAC_MAX_ENTRY}; "
                          "using the closed-form bound", stacklevel=2)
        else:
            try:
                record = solve_vac(merged, timeout=timeout)
            except (OSError, VacError) as error:
                warnings.warn(f"vac is unavailable ({error}); using the "
                              "closed-form bound", stacklevel=2)
            else:
                proven = record["proven_optimal"]
                vac_bound = record["length"] if proven else record["lower_bound"]
                # Merging units loses exactness, so only unmerged runs are proven.
                info.update(source="vac", vac=record,
                            proven_optimal=proven and merged_units is units)
                bound = max(bound, vac_bound)

    return (bound, info) if return_info else bound
