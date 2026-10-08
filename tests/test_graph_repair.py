import random

import networkx as nx
import pytest
from networkx.algorithms.isomorphism import categorical_edge_match, categorical_node_match

import assemblycfg as cfg
from test_det import AMINO_ACIDS

# Lipids from examples/lipids.py with their assembly indices.
LIPIDS = (
    ("propanoic-acid", "CCC(=O)O", 3),
    ("nonanoic-acid", "CCCCCCCCC(=O)O", 5),
    ("triglyceride-2", "CCC(=O)OCC(COC(=O)CC)OC(=O)CC", 7),
    ("phospholipid-2", "CCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CC", 14),
)


def fragment(graph, result, edges):
    return graph.edge_subgraph(result.edges[i] for i in edges)


def isomorphic(first, second):
    return nx.is_isomorphic(first, second,
                            node_match=categorical_node_match('color', None),
                            edge_match=categorical_edge_match('color', None))


def assert_valid_certificate(graph, result):
    """Replay the construction: every join is real and the residual tiles the graph."""
    assert len(result.edges) == graph.number_of_edges()
    assert set(map(frozenset, result.edges)) == set(map(frozenset, graph.edges))
    symbols = {}
    for terminal in result.terminals:
        assert len(terminal.edges) == 1
        assert terminal.symbol not in symbols
        symbols[terminal.symbol] = fragment(graph, result, terminal.edges)
    for rule in result.rules:
        assert rule.left in symbols and rule.right in symbols and rule.id not in symbols
        assert not set(rule.left_edges) & set(rule.right_edges)
        assert set(rule.left_edges) | set(rule.right_edges) == set(rule.edges)
        left = fragment(graph, result, rule.left_edges)
        right = fragment(graph, result, rule.right_edges)
        assert set(left) & set(right), "nonincident join"
        assert isomorphic(left, symbols[rule.left])
        assert isomorphic(right, symbols[rule.right])
        symbols[rule.id] = fragment(graph, result, rule.edges)
    covered = [edge for residual in result.residual for edge in residual.edges]
    assert sorted(covered) == list(range(len(result.edges)))
    for residual in result.residual:
        piece = fragment(graph, result, residual.edges)
        assert nx.is_connected(piece)
        assert isomorphic(piece, symbols[residual.symbol])
    adjustment = result.components if result.compensate_disjoint else 1
    assert result.upper_bound == result.rule_count + result.remaining_fragments - adjustment
    assert result.upper_bound <= result.trivial_upper_bound


@pytest.mark.parametrize("compensate_disjoint", [False, True])
@pytest.mark.parametrize(
    "smiles",
    [smiles for _, smiles, _ in LIPIDS + AMINO_ACIDS]
    + ["c1ccc2ccccc2c1", "C1CC2CCC1CC2", "CC(C)(C)C(C)(C)C", "CCO.CCO.O"],
)
def test_certificates_replay(smiles, compensate_disjoint):
    graph = cfg.smi_to_nx(smiles, add_hydrogens=False)
    assert_valid_certificate(graph, cfg.graph_repair(graph, compensate_disjoint=compensate_disjoint))


@pytest.mark.parametrize(("smiles", "assembly_index"), [(s, a) for _, s, a in LIPIDS])
def test_lipids_reach_their_assembly_index(smiles, assembly_index):
    graph = cfg.smi_to_nx(smiles, add_hydrogens=False)
    assert cfg.calculate_assembly_path_graph_repair(graph)[0] == assembly_index


@pytest.mark.parametrize(("smiles", "lower_bound"), [(s, a) for _, s, a in AMINO_ACIDS])
def test_amino_acid_upper_bounds(smiles, lower_bound):
    graph = cfg.smi_to_nx(smiles, add_hydrogens=False)
    result = cfg.graph_repair(graph)
    assert lower_bound <= result.upper_bound <= result.trivial_upper_bound


@pytest.mark.parametrize(("smiles", "expected"), [("CCO", 1), ("CCCC", 2), ("CCCCCC", 3), ("c1ccccc1", 3)])
def test_small_molecules(smiles, expected):
    assert cfg.graph_repair(cfg.smi_to_nx(smiles, add_hydrogens=False)).upper_bound == expected


def test_reuses_branched_fragments():
    # Two tert-butyl groups: the C(C)(C)C star is built once and reused.
    graph = cfg.smi_to_nx("CC(C)(C)C(C)(C)C", add_hydrogens=False)
    assert cfg.graph_repair(graph).upper_bound == 4


def test_disconnected_conventions():
    graph = cfg.smi_to_nx("CCO.CCO", add_hydrogens=False)
    default = cfg.graph_repair(graph)
    compensated = cfg.graph_repair(graph, compensate_disjoint=True)
    assert default.components == compensated.components == 2
    assert default.upper_bound == compensated.upper_bound + 1 == 2


def test_graphs_without_bonds():
    graph = nx.Graph()
    assert cfg.calculate_assembly_path_graph_repair(graph) == (0, None, None)
    graph.add_node(0, color="C")
    assert cfg.graph_repair(graph).upper_bound == 0
    graph.add_edge(1, 2, color=1)
    graph.nodes[1]["color"] = graph.nodes[2]["color"] = "O"
    result = cfg.graph_repair(graph)
    assert (result.upper_bound, result.components) == (0, 1)


def test_pathway_outputs():
    graph = cfg.smi_to_nx(LIPIDS[2][1], add_hydrogens=False)
    length, virtual_objects, path = cfg.calculate_assembly_path_graph_repair(graph)
    result = cfg.graph_repair(graph)
    assert length == result.upper_bound
    assert len(virtual_objects) == path.number_of_nodes() == len(result.terminals) + result.rule_count
    assert nx.is_directed_acyclic_graph(path)
    for rule in result.rules:
        assert set(path.successors(rule.id)) == {rule.left, rule.right}
        assert virtual_objects[rule.id].number_of_edges() == len(rule.edges)
    for virtual_object in virtual_objects:
        assert set(virtual_object.edges) <= set(graph.edges)


def test_deterministic():
    graph = cfg.smi_to_nx(LIPIDS[3][1], add_hydrogens=False)
    assert cfg.graph_repair(graph) == cfg.graph_repair(graph)
    assert cfg.graph_repair(graph, iterations=5, seed=1) == cfg.graph_repair(graph, iterations=5, seed=1)


@pytest.mark.parametrize("iterations", [0, -1])
def test_rejects_non_positive_iterations(iterations):
    with pytest.raises(ValueError, match="iterations"):
        cfg.graph_repair(nx.Graph(), iterations=iterations)


def test_more_iterations_do_not_worsen_the_bound():
    graph = nx.disjoint_union_all(
        [cfg.smi_to_nx(smiles, add_hydrogens=False) for _, smiles, _ in AMINO_ACIDS]
    )
    single = cfg.graph_repair(graph)
    repeated = cfg.graph_repair(graph, iterations=10, seed=0)
    assert repeated.upper_bound <= single.upper_bound
    assert_valid_certificate(graph, repeated)


@pytest.mark.parametrize("seed", range(5))
def test_bound_survives_shuffled_atom_order(seed):
    # Counting pairs greedily in bond order misses repeats along shuffled
    # chains; maximum matching does not.
    graph = cfg.smi_to_nx("C" * 16 + "C(=O)O", add_hydrogens=False)
    rng = random.Random(seed)
    nodes, edges = list(graph.nodes), list(graph.edges)
    rng.shuffle(nodes)
    rng.shuffle(edges)
    shuffled = nx.Graph()
    shuffled.add_nodes_from((node, graph.nodes[node]) for node in nodes)
    shuffled.add_edges_from((u, v, graph.edges[u, v]) for u, v in edges)
    result = cfg.graph_repair(shuffled)
    assert_valid_certificate(shuffled, result)
    assert result.upper_bound == 6
