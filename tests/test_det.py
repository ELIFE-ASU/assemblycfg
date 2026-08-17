import random
from collections import Counter

import networkx as nx
import pytest

import assemblycfg as cfg
import assemblycfg.det as det

TRYPTOPHAN = "c1[nH]c2ccccc2c1C[C@H](N)C(=O)O"

AMINO_ACIDS = (
    ("glycine", "C(C(=O)O)N", 3),
    ("alanine", "C[C@@H](C(=O)O)N", 4),
    ("serine", "C([C@@H](C(=O)O)N)O", 4),
    ("aspartic-acid", "C([C@@H](C(=O)O)N)C(=O)O", 5),
    ("cysteine", "C([C@@H](C(=O)O)N)S", 5),
    ("glutamic-acid", "C(CC(=O)O)[C@@H](C(=O)O)N", 5),
    ("threonine", "C[C@H]([C@@H](C(=O)O)N)O", 5),
    ("valine", "CC(C)[C@@H](C(=O)O)N", 5),
    ("asparagine", "C([C@@H](C(=O)O)N)C(=O)N", 6),
    ("glutamine", "C(CC(=O)N)[C@@H](C(=O)O)N", 6),
    ("isoleucine", "CC[C@H](C)[C@@H](C(=O)O)N", 6),
    ("leucine", "CC(C)C[C@@H](C(=O)O)N", 6),
    ("lysine", "C(CCN)C[C@@H](C(=O)O)N", 6),
    ("proline", "C1C[C@H](NC1)C(=O)O", 6),
    ("methionine", "CSCC[C@@H](C(=O)O)N", 7),
    ("arginine", "C(C[C@@H](C(=O)O)N)CN=C(N)N", 9),
    ("histidine", "C1=C(NC=N1)C[C@@H](C(=O)O)N", 9),
    ("phenylalanine", "C1=CC=C(C=C1)C[C@@H](C(=O)O)N", 9),
    ("tyrosine", "C1=CC(=CC=C1C[C@@H](C(=O)O)N)O", 9),
    ("tryptophan", TRYPTOPHAN, 11),
)

# Multiplicity is retained because this fixture represents a joint target, not a
# set of unique SMILES strings.
ATMOSPHERIC_NETWORK_SMILES = (
    "O=O",
    "O",
    "[OH]",
    "OO[H]",
    "OO",
    "[HH]",
    "[C]=O",
    "[CH]=O",
    "C=O",
    "C",
    "[CH3]",
    "CC",
    "[N]=O",
    "N(=O)[O]",
    "[NH]=O",
    "[O-][O+]=O",
    "O=[N](=O)O",
    "C=C=C",
    "[CH]=C=C",
    "CC#C",
    "C=C=C",
    "C=CC[CH2]",
    "CCC=O",
    "CC=C",
    "CCC[CH2]",
    "CCC",
    "CCO",
    "C=CO",
    "CC[CH2]",
    "C=C",
    "[CH]",
    "COO[CH2]",
    "C[O]",
    "C=C=O",
    "CC(=O)[CH2]",
    "CC=O",
    "C#C",
    "[CH]=C",
    "[C]#[C]",
    "[CH]=C",
    "[CH]=S",
    "S=C=S",
    "[C]=S",
    "O=C=S",
    "[SH]",
    "S",
    "O=S(=O)=O",
    "S=S",
    "O=[SH]",
    "O=S(=O)(O)O",
    "O=S=O",
    "[S]=O",
    "O=C=O",
    "N#N",
)


@pytest.fixture
def seed_det(monkeypatch):
    def apply_seed(seed=0):
        rng = random.Random(seed)
        monkeypatch.setattr(det.random, "choice", rng.choice)

    return apply_seed


def assert_valid_det_result(graph, result, lower_bound, upper_bound):
    path_length, virtual_objects, path = result

    assert lower_bound <= path_length <= upper_bound
    assert virtual_objects is not None
    assert path is not None
    assert len(virtual_objects) == path.number_of_nodes()
    assert set(path) == set(range(len(virtual_objects)))
    assert nx.is_directed_acyclic_graph(path)
    assert not list(nx.isolates(path))

    for virtual_object in virtual_objects:
        assert isinstance(virtual_object, nx.Graph)
        assert virtual_object.number_of_edges() >= 1


def test_partition_into_disjoint_trails_covers_every_edge_once(seed_det):
    graph = nx.Graph([(0, 1), (1, 2), (2, 0), (1, 3), (3, 4)])
    seed_det()

    trails = det.partition_into_disjoint_trails(graph)
    traversed_edges = Counter(
        frozenset((left, right))
        for trail in trails
        for left, right in zip(trail, trail[1:])
    )

    assert traversed_edges == Counter(frozenset(edge) for edge in graph.edges)
    assert sum(len(trail) - 1 for trail in trails) == graph.number_of_edges()


def test_partition_into_disjoint_trails_handles_empty_graph():
    assert det.partition_into_disjoint_trails(nx.Graph()) == []


def test_trails_to_sequences_encodes_units_and_counts_single_edges():
    graph = nx.Graph()
    graph.add_nodes_from(
        [(0, {"color": "C"}), (1, {"color": "O"}), (2, {"color": "N"})]
    )
    graph.add_nodes_from([(3, {"color": "S"}), (4, {"color": "H"})])
    graph.add_edge(0, 1, color=1)
    graph.add_edge(1, 2, color=2)
    graph.add_edge(3, 4, color=1)

    sequences, unit_to_symbol, trivial_count = det.trails_to_sequences(
        [[0, 1, 2], [3, 4]], graph
    )

    assert sequences == ["\x01\x02"]
    assert unit_to_symbol == {("C", 1, "O"): "\x01", ("O", 2, "N"): "\x02"}
    assert trivial_count == 1


def test_trails_to_sequences_rejects_a_missing_edge():
    graph = nx.Graph()
    graph.add_nodes_from([(0, {"color": "C"}), (1, {"color": "O"})])

    with pytest.raises(KeyError, match=r"Edge \(0, 1\) not found"):
        det.trails_to_sequences([[0, 1, 0]], graph)


@pytest.mark.parametrize(
    ("sequences", "expected_sequences", "expected_rules"),
    [
        pytest.param([], [], [], id="empty"),
        pytest.param([""], [[]], [], id="empty-sequence"),
        pytest.param(
            ["aaaa"],
            [["NT_0", "NT_0"]],
            [("NT_0", ("a", "a"))],
            id="overlapping-pair",
        ),
        pytest.param(
            ["ababcdcd"],
            [["NT_0", "NT_0", "NT_1", "NT_1"]],
            [("NT_0", ("a", "b")), ("NT_1", ("c", "d"))],
            id="two-pairs",
        ),
    ],
)
def test_repair_compression(sequences, expected_sequences, expected_rules):
    assert det.repair_compression(sequences) == (expected_sequences, expected_rules)


def test_repair_compression_rules_round_trip():
    original = "ababcdcd"
    final_sequences, rules = det.repair_compression([original])

    rebuilt = "".join(det.unpack_path(symbol, rules) for symbol in final_sequences[0])

    assert rebuilt == original


def test_compute_assembly_path_length_includes_each_term():
    assert (
        det.compute_assembly_path_length_from_compression(
            [["NT_0", "a"], ["b"]], [("NT_0", ("a", "a"))], ai_correction=2
        )
        == 6
    )


def test_purge_unique_units_keeps_repeated_reversible_units():
    graph = nx.Graph()
    graph.add_nodes_from(
        [
            (0, {"color": "C"}),
            (1, {"color": "O"}),
            (2, {"color": "O"}),
            (3, {"color": "C"}),
            (4, {"color": "N"}),
            (5, {"color": "N"}),
        ]
    )
    graph.add_edges_from([(0, 1, {"color": 1}), (2, 3, {"color": 1})])
    graph.add_edge(4, 5, color=2)

    purged, unique_count = det.purge_unique_units(graph)

    assert unique_count == 1
    assert set(purged) == {0, 1, 2, 3}
    assert set(purged.edges) == {(0, 1), (2, 3)}
    assert graph.number_of_edges() == 3


def test_process_paths_uses_one_node_per_symbol():
    rules = [
        ("NT_0", ("a", "b")),
        ("NT_1", ("NT_0", "a")),
    ]
    unit_to_symbol = {("C", 1, "O"): "a", ("O", 2, "N"): "b"}

    virtual_objects, path = det.process_paths(rules, unit_to_symbol)

    assert set(path) == {0, 1, 2, 3}
    assert set(path.edges) == {(0, 1), (0, 2), (3, 0), (3, 1)}
    assert [graph.number_of_edges() for graph in virtual_objects] == [2, 1, 1, 3]
    assert not list(nx.isolates(path))


def test_det_returns_zero_for_trivial_graphs_without_output(capsys):
    one_edge = nx.Graph()
    one_edge.add_edge(0, 1, color=1)
    nx.set_node_attributes(one_edge, "C", "color")

    disconnected_edges = nx.disjoint_union(one_edge, one_edge)

    assert cfg.calculate_assembly_path_det(nx.Graph()) == (0, None, None)
    assert cfg.calculate_assembly_path_det(one_edge) == (0, None, None)
    assert cfg.calculate_assembly_path_det(disconnected_edges) == (0, None, None)
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("iterations", [0, -1])
def test_det_rejects_non_positive_iterations(iterations):
    with pytest.raises(ValueError, match="iterations must be at least 1"):
        cfg.calculate_assembly_path_det(nx.Graph(), iterations=iterations)


@pytest.mark.parametrize(
    ("smiles", "add_hydrogens", "lower_bound", "accepted_upper_bound"),
    [
        pytest.param("C1CCC1", True, 4, 9, id="cyclobutane-with-hydrogen"),
        pytest.param("C1CCC1", False, 2, 3, id="cyclobutane"),
        pytest.param("C" * 52, False, 7, 9, id="long-alkane"),
        pytest.param("c1ccccc1", True, 5, 11, id="benzene"),
        pytest.param(TRYPTOPHAN, False, 11, 14, id="tryptophan"),
    ],
)
def test_known_molecule_upper_bounds(
    smiles, add_hydrogens, lower_bound, accepted_upper_bound, seed_det
):
    graph = cfg.smi_to_nx(smiles, add_hydrogens=add_hydrogens)
    seed_det()

    result = cfg.calculate_assembly_path_det(graph)

    assert_valid_det_result(graph, result, lower_bound, accepted_upper_bound)


def test_small_joint_system_upper_bound(seed_det):
    graphs = [cfg.smi_to_nx(smiles) for smiles in ("O=C=O", "C", "O", "N")]
    graph = nx.disjoint_union_all(graphs)
    seed_det()

    result = cfg.calculate_assembly_path_det(graph)

    assert_valid_det_result(graph, result, lower_bound=6, upper_bound=7)


def test_more_iterations_do_not_worsen_the_best_path(seed_det):
    graph = cfg.smi_to_nx(TRYPTOPHAN, add_hydrogens=False)
    seed_det()
    single_iteration = cfg.calculate_assembly_path_det(graph, iterations=1)
    seed_det()
    repeated_iterations = cfg.calculate_assembly_path_det(graph, iterations=100)

    assert repeated_iterations[0] <= single_iteration[0]
    assert_valid_det_result(
        graph,
        repeated_iterations,
        lower_bound=11,
        upper_bound=single_iteration[0],
    )


@pytest.mark.parametrize(
    ("smiles", "known_lower_bound"),
    [(smiles, lower_bound) for _, smiles, lower_bound in AMINO_ACIDS],
    ids=[name for name, _, _ in AMINO_ACIDS],
)
def test_amino_acid_upper_bounds(smiles, known_lower_bound, seed_det):
    graph = cfg.smi_to_nx(smiles, add_hydrogens=False)
    seed_det()

    result = cfg.calculate_assembly_path_det(graph)

    assert_valid_det_result(
        graph,
        result,
        lower_bound=known_lower_bound,
        upper_bound=graph.number_of_edges(),
    )


def test_joint_amino_acid_upper_bound(seed_det):
    graphs = [
        cfg.smi_to_nx(smiles, add_hydrogens=False) for _, smiles, _ in AMINO_ACIDS
    ]
    joint_graph = nx.disjoint_union_all(graphs)
    seed_det()

    result = cfg.calculate_assembly_path_det(joint_graph, iterations=100)

    assert_valid_det_result(
        joint_graph,
        result,
        lower_bound=37,
        upper_bound=joint_graph.number_of_edges(),
    )


def test_large_disconnected_joint_system(seed_det):
    graphs = [cfg.smi_to_nx(smiles) for smiles in ATMOSPHERIC_NETWORK_SMILES]
    joint_graph = nx.disjoint_union_all(graphs)
    seed_det()

    result = cfg.calculate_assembly_path_det(joint_graph)

    assert_valid_det_result(
        joint_graph,
        result,
        lower_bound=3,
        upper_bound=joint_graph.number_of_edges(),
    )
