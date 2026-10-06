import random
from collections import Counter

import networkx as nx
import pytest

import assemblycfg as cfg
from assemblycfg.cfg_ai import convert_to_cnf, repair, repair_pairs

ABRACADABRA_OBJECTS = {
    "a",
    "b",
    "c",
    "d",
    "r",
    "ab",
    "abr",
    "abra",
    "abrac",
    "abraca",
    "abracad",
}


@pytest.mark.parametrize(
    ("source", "expected_index"),
    [
        pytest.param("abracadabra", 7, id="abracadabra"),
        pytest.param("aaaaaaa", 4, id="repeated-symbol"),
        pytest.param("ababcdcd", 5, id="two-repeated-pairs"),
        pytest.param("aaaadbbbbcaa", 8, id="two-symbol-runs"),
        pytest.param(["aaaa", "bbbbcaa"], 6, id="joint-input"),
    ],
)
def test_known_assembly_indices(source, expected_index):
    assembly_index, _ = cfg.ai_core(source)

    assert assembly_index == expected_index


def test_repair_shares_productions_across_joint_input():
    symbols, productions = repair(["aaaa", "bbbbcaa"])

    assert symbols == [["A1", "A1"], ["A2", "A2", "c", "A1"]]
    assert productions == {"A1": ["a", "a"], "A2": ["b", "b"]}



def naive_repair(sequences):
    """Recount every pair each round: the quadratic definition of RePair."""
    seqs = [list(seq) for seq in sequences]
    rules = []
    while True:
        counts = Counter()
        for seq in seqs:
            counts.update(zip(seq, seq[1:]))
        frequent = {pair: c for pair, c in counts.items() if c > 1}
        if not frequent:
            return seqs, rules
        pair = max(frequent, key=frequent.get)
        rules.append((f"R{len(rules)}", pair))
        for idx, seq in enumerate(seqs):
            out, i = [], 0
            while i < len(seq):
                if tuple(seq[i:i + 2]) == pair:
                    out.append(rules[-1][0])
                    i += 2
                else:
                    out.append(seq[i])
                    i += 1
            seqs[idx] = out


@pytest.mark.parametrize("alphabet", ["a", "ab", "abc", "acgt"])
def test_repair_pairs_matches_naive_repair(alphabet):
    rng = random.Random(alphabet)
    for _ in range(200):
        sequences = ["".join(rng.choice(alphabet) for _ in range(rng.randint(0, 40)))
                     for _ in range(rng.randint(1, 3))]
        assert repair_pairs(sequences, lambda i: f"R{i}") == naive_repair(sequences)

@pytest.mark.parametrize("source", ["Uppercase", "lowercase1"])
def test_repair_rejects_non_lowercase_ascii(source):
    with pytest.raises(
        AssertionError,
        match="Input string must consist of lowercase ASCII characters only",
    ):
        repair(source)


def test_repair_rejects_non_string_joint_member():
    with pytest.raises(TypeError, match="Input must be a string or a list of strings"):
        repair(["valid", 1])


def test_convert_to_cnf_splits_long_productions_and_start_symbol():
    start, productions = convert_to_cnf("abc", {"A": ["a", "b", "c"]})

    assert start == "S_0"
    assert productions["T_a"] == ["a"]
    assert productions["T_b"] == ["b"]
    assert productions["T_c"] == ["c"]
    assert productions["N1"] == ["T_a", "T_b"]
    assert productions["A"] == ["N1", "T_c"]
    assert productions["N2"] == ["T_a", "T_b"]
    assert productions["S_0"] == ["N2", "T_c"]
    assert all(1 <= len(expansion) <= 2 for expansion in productions.values())


def test_abracadabra_pathway_has_expected_objects_and_dependencies():
    assembly_index, virtual_objects, path = cfg.repair_with_pathways("abracadabra")

    assert assembly_index == 7
    assert Counter(virtual_objects) == Counter(ABRACADABRA_OBJECTS)
    assert set(path) == ABRACADABRA_OBJECTS | {"abracadabra"}
    assert set(path.successors("abracad")) == {"abracadabra"}
    assert set(path.successors("abra")) == {"abrac", "abracadabra"}
    assert {node for node, degree in path.out_degree if degree == 0} == {"abracadabra"}
    assert nx.is_directed_acyclic_graph(path)

    for component, product in path.edges:
        assert component in product
        assert len(component) < len(product)


def test_joint_pathway_contains_both_products():
    assembly_index, virtual_objects, path = cfg.repair_with_pathways(
        ["aaaa", "bbbbcaa"]
    )

    assert assembly_index == 6
    assert virtual_objects
    assert {"aaaa", "bbbbcaa"} <= set(path)
    assert {node for node, degree in path.out_degree if degree == 0} == {
        "aaaa",
        "bbbbcaa",
    }
    assert nx.is_directed_acyclic_graph(path)
