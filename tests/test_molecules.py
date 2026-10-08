from collections import Counter

import networkx as nx
import pytest
from rdkit import Chem

import assemblycfg as cfg


@pytest.fixture
def ethanol_molfile(tmp_path):
    path = tmp_path / "ethanol.mol"
    Chem.MolToMolFile(Chem.MolFromSmiles("CCO"), str(path))
    return path


@pytest.mark.parametrize(
    ("bond_type", "expected_order"),
    [
        (Chem.BondType.SINGLE, 1),
        (Chem.BondType.DOUBLE, 2),
        (Chem.BondType.TRIPLE, 3),
        (Chem.BondType.QUADRUPLE, 4),
        (Chem.BondType.QUINTUPLE, 5),
        (Chem.BondType.IONIC, 6),
        pytest.param(Chem.BondType.AROMATIC, 1, id="unknown-default"),
    ],
)
def test_bond_order_conversion(bond_type, expected_order):
    assert cfg.bond_order_rdkit_to_int(bond_type) == expected_order


def test_smiles_to_graph_controls_explicit_hydrogens():
    with_hydrogens = cfg.smi_to_nx("C")
    without_hydrogens = cfg.smi_to_nx("C", add_hydrogens=False)

    assert Counter(nx.get_node_attributes(with_hydrogens, "color").values()) == {
        "C": 1,
        "H": 4,
    }
    assert with_hydrogens.number_of_edges() == 4
    assert nx.get_node_attributes(without_hydrogens, "color") == {0: "C"}
    assert without_hydrogens.number_of_edges() == 0


def test_smiles_to_graph_preserves_bond_order():
    graph = cfg.smi_to_nx("C#N", add_hydrogens=False)

    assert Counter(nx.get_node_attributes(graph, "color").values()) == {"C": 1, "N": 1}
    assert list(nx.get_edge_attributes(graph, "color").values()) == [3]


def test_disconnected_smiles_emits_a_warning():
    with pytest.warns(UserWarning, match="Disconnected molecules detected"):
        molecule = cfg.smi_to_mol("C.O", add_hydrogens=False)

    assert len(Chem.GetMolFrags(molecule)) == 2


def test_smiles_can_be_parsed_without_sanitization():
    molecule = cfg.smi_to_mol("C", sanitize=False)

    assert isinstance(molecule, Chem.Mol)
    assert molecule.GetNumAtoms() == 1


@pytest.mark.parametrize("converter", [cfg.smi_to_mol, cfg.smi_to_nx])
def test_smiles_converters_reject_invalid_input(converter):
    with pytest.raises(ValueError, match="Unable to parse SMILES"):
        converter("not a SMILES")


def test_graph_dictionary_round_trip_preserves_topology_and_attributes():
    graph = nx.Graph()
    graph.add_node("carbon", color="C")
    graph.add_node("oxygen", color="O")
    graph.add_node("defaulted")
    graph.add_edge("carbon", "oxygen", color=2)
    graph.add_edge("oxygen", "defaulted")

    graph_dict = cfg.nx_to_dict(graph)
    rebuilt = cfg.dict_to_nx(graph_dict)

    assert graph_dict["defaulted"]["vertex_color"] == "C"
    assert graph_dict["oxygen"]["edge_colors"] == [2, 1]
    assert nx.get_node_attributes(rebuilt, "color") == {
        "carbon": "C",
        "oxygen": "O",
        "defaulted": "C",
    }
    assert nx.get_edge_attributes(rebuilt, "color") == {
        ("carbon", "oxygen"): 2,
        ("oxygen", "defaulted"): 1,
    }


@pytest.mark.parametrize(
    ("graph_dict", "message"),
    [
        pytest.param(
            {0: {"vertex_color": "C", "neighbors": [1], "edge_colors": []}},
            "different numbers of neighbors and edge colors",
            id="misaligned-edge-colors",
        ),
        pytest.param(
            {0: {"vertex_color": "C", "neighbors": [1], "edge_colors": [1]}},
            "references unknown neighbor",
            id="unknown-neighbor",
        ),
        pytest.param(
            {
                0: {"vertex_color": "C", "neighbors": [1], "edge_colors": [1]},
                1: {"vertex_color": "O", "neighbors": [0], "edge_colors": [2]},
            },
            "has conflicting colors",
            id="conflicting-reciprocal-colors",
        ),
    ],
)
def test_dict_to_nx_rejects_inconsistent_adjacency(graph_dict, message):
    with pytest.raises(ValueError, match=message):
        cfg.dict_to_nx(graph_dict)


def test_remove_hydrogen_mutates_and_returns_the_input_graph():
    graph = nx.Graph()
    graph.add_nodes_from([(0, {"color": "C"}), (1, {"color": "H"})])
    graph.add_edge(0, 1, color=1)

    result = cfg.remove_hydrogen_from_graph(graph)

    assert result is graph
    assert list(graph) == [0]
    assert graph.number_of_edges() == 0


def test_get_disconnected_subgraphs_returns_each_component():
    graph = nx.Graph([(0, 1), (2, 3)])

    components = cfg.get_disconnected_subgraphs(graph)

    assert {frozenset(component) for component in components} == {
        frozenset({0, 1}),
        frozenset({2, 3}),
    }


@pytest.mark.parametrize("safe_sanitise", [False, True], ids=["strict", "safe"])
def test_molfile_to_mol_supports_both_sanitizers(ethanol_molfile, safe_sanitise):
    molecule = cfg.molfile_to_mol(
        str(ethanol_molfile),
        add_hydrogens=False,
        safe_sanitise=safe_sanitise,
    )

    assert isinstance(molecule, Chem.Mol)
    assert Counter(atom.GetSymbol() for atom in molecule.GetAtoms()) == {"C": 2, "O": 1}


def test_mol2graph_returns_a_hydrogen_free_dictionary(ethanol_molfile):
    graph_dict = cfg.mol2graph(str(ethanol_molfile))

    assert Counter(node["vertex_color"] for node in graph_dict.values()) == {
        "C": 2,
        "O": 1,
    }
    assert sum(len(node["neighbors"]) for node in graph_dict.values()) == 4
    assert {
        edge_color for node in graph_dict.values() for edge_color in node["edge_colors"]
    } == {1}


def test_print_helpers_render_graph_dictionaries(capsys):
    graph = nx.Graph()
    graph.add_node(0, color="C")

    cfg.print_virtual_objects([graph])

    assert capsys.readouterr().out == (
        "0: \n0: {'vertex_color': 'C', 'neighbors': [], 'edge_colors': []}\n\n"
    )
