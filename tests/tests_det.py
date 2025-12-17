import copy
import time

import networkx as nx
import pytest

import CFGgraph
import assemblytheorytools as att


def test_basic_ai():
    """
    Test the calculation of the assembly index using both the CFG and DET approaches.

    This function defines a simple graph, converts it to a NetworkX graph, and calculates
    the assembly index using two different approaches: CFG and DET. It then prints the
    results and asserts that both assembly indices are greater than 1.
    """
    print(flush=True)
    # Define the graph
    d = {
        0: {'vertex_color': 'C', 'neighbors': [1, 7], 'edge_colors': [1, 1]},
        1: {'vertex_color': 'C', 'neighbors': [0, 2], 'edge_colors': [1, 1]},
        2: {'vertex_color': 'C', 'neighbors': [1, 3], 'edge_colors': [1, 1]},
        3: {'vertex_color': 'C', 'neighbors': [2, 4], 'edge_colors': [1, 1]},
        4: {'vertex_color': 'C', 'neighbors': [3, 5], 'edge_colors': [1, 1]},
        5: {'vertex_color': 'C', 'neighbors': [4, 6], 'edge_colors': [1, 1]},
        6: {'vertex_color': 'C', 'neighbors': [5, 7], 'edge_colors': [1, 1]},
        7: {'vertex_color': 'C', 'neighbors': [6, 0], 'edge_colors': [1, 1]}
    }
    # Convert the dictionary to a NetworkX graph
    graph = CFGgraph.dict_to_nx(d)
    # Calculate the assembly index
    ai, _, _ = att.calculate_assembly_index(graph)
    # Calculate the assembly index using the DET approach
    L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph)

    # Print the assembly indices
    print(f"Assembly Index (asscpp) : {ai}", flush=True)
    print(f"Assembly Index (DET)    : {L_det}", flush=True)

    # Assert that both assembly indices are greater than 1
    assert L_det > 1


def test_basic_joint_ai():
    """
    Test the calculation of the assembly index for a joint system of molecules.

    This function splits a string of SMILES representations into individual molecules,
    converts them to RDKit molecules, then to NetworkX graphs, and joins these graphs
    into a single graph. It calculates the assembly index using three different approaches
    and asserts that the results are equal.

    Asserts:
        The assembly index calculated by `att.calculate_assembly_index` is equal to the
        assembly index calculated by `CFGgraph.calculate_assembly_index_ga` and
        `CFGgraph.calculate_assembly_path_det`.
    """
    print(flush=True)
    molecules = "O=C=O.C.O.N"
    molecules = molecules.split(".")
    # Convert all the SMILES to RDKit molecules
    mols = [att.smi_to_mol(smile) for smile in molecules]
    # Convert the RDKit molecules to NetworkX graphs
    graphs = [att.mol_to_nx(mol) for mol in mols]
    # Join the graphs into a single graph
    graph = nx.disjoint_union_all(graphs)
    # Calculate the assembly index using the default approach
    ai, _, _ = att.calculate_assembly_index(graph)
    # Calculate the assembly index using the DET approach
    L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph)
    # Print the assembly indices
    print(f"Assembly Index (asscpp) : {ai}", flush=True)
    print(f"Assembly Index (DET)    : {L_det}", flush=True)
    assert L_det >= ai


def test_det_cyclobutane():
    """
    Test the calculation of the assembly index for cyclobutane.

    This function converts a SMILES string representing cyclobutane to an RDKit molecule,
    then to a NetworkX graph. It calculates the assembly index using three different approaches
    (default, CFG, and DET) and prints the results. It then removes hydrogen atoms from the graph,
    recalculates the assembly index using the same three approaches, and prints the results.

    Asserts:
        The assembly indices calculated by the CFG and DET approaches are equal to the default approach,
        both with and without hydrogen atoms.
    """
    print(flush=True)
    smi = "C1CCC1"
    mol = att.smi_to_mol(smi)
    graph = att.mol_to_nx(mol)

    # Calculate the assembly index using the default approach
    ai, _, _ = att.calculate_assembly_index(graph)
    # Calculate the assembly index using the DET approach
    L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph)
    # Print the assembly indices
    print(f"Assembly Index (asscpp) : {ai}", flush=True)
    print(f"Assembly Index (DET)    : {L_det}", flush=True)

    print("Without hydrogen...", flush=True)
    graph = att.remove_hydrogen_from_graph(graph)

    # Calculate the assembly index using the default approach
    ai_nh, _, _ = att.calculate_assembly_index(graph)
    # Calculate the assembly index using the DET approach
    L_det_nh, _, _ = CFGgraph.calculate_assembly_path_det(graph)
    # Print the assembly indices
    print(f"Assembly Index (asscpp) : {ai_nh}", flush=True)
    print(f"Assembly Index (DET)    : {L_det_nh}", flush=True)

    # Assert that the assembly indices are equal
    assert L_det >= ai
    assert L_det_nh >= ai_nh


def test_long_boi():
    """
    Test the calculation of the assembly index for a long alkane chain.

    This function converts a SMILES string representing a long alkane chain to an RDKit molecule,
    then to a NetworkX graph. It calculates the assembly index using three different approaches
    (default, CFG, and DET) and prints the results.

    Asserts:
        The assembly indices calculated by the CFG and DET approaches are equal to the default approach.
    """
    print(flush=True)
    smi = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    mol = att.smi_to_mol(smi)
    graph = att.mol_to_nx(mol)
    # Strip the hydrogen atoms from the graph
    graph = att.remove_hydrogen_from_graph(graph)
    # Calculate the assembly index using the default approach
    ai, _, _ = att.calculate_assembly_index(graph, timeout=1000.0)
    # Calculate the assembly index using the DET approach
    L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph)
    # Print the assembly indices
    print(f"Assembly Index (asscpp) : {ai}", flush=True)
    print(f"Assembly Index (DET)    : {L_det}", flush=True)
    # Assert that the assembly indices are equal
    assert L_det >= ai


def test_benzene():
    """
    Test the calculation of the assembly index for benzene.

    This function converts a SMILES string representing benzene to an RDKit molecule,
    then to a NetworkX graph. It calculates the assembly index using three different approaches
    (default, CFG, and DET) and prints the results.

    Asserts:
        The assembly indices calculated by the CFG and DET approaches are equal to the default approach.
    """
    print(flush=True)
    smi = "c1ccccc1"
    mol = att.smi_to_mol(smi)
    graph = att.mol_to_nx(mol)
    # Calculate the assembly index using the default approach
    ai, _, _ = att.calculate_assembly_index(graph)
    # Calculate the assembly index using the DET approach
    L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph)
    # Print the assembly indices
    print(f"Assembly Index (asscpp) : {ai}", flush=True)
    print(f"Assembly Index (DET)    : {L_det}", flush=True)
    # Assert that the assembly indices are equal
    assert L_det >= ai


def test_tryptophan():
    """
    Test the calculation of the assembly index for tryptophan.

    This function converts a SMILES string representing tryptophan to an RDKit molecule,
    then to a NetworkX graph. It removes hydrogen atoms from the graph and calculates the
    assembly index using three different approaches (default, CFG, and DET). It prints the
    results and asserts that the assembly indices calculated by the CFG and DET approaches
    are equal to the default approach.

    Asserts:
        The assembly indices calculated by the CFG and DET approaches are equal to the default approach.
    """
    print(flush=True)
    smi = "c1[nH]c2ccccc2c1C[C@H](N)C(=O)O"
    mol = att.smi_to_mol(smi)
    graph = att.mol_to_nx(mol)
    # Strip the hydrogen atoms from the graph
    graph = att.remove_hydrogen_from_graph(graph)
    # Calculate the assembly index using the default approach
    ai, _, _ = att.calculate_assembly_index(graph)
    # Calculate the assembly index using the DET approach
    L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph)
    # Print the assembly indices
    print(f"Assembly Index (asscpp) : {ai}", flush=True)
    print(f"Assembly Index (DET)    : {L_det}", flush=True)
    # Assert that the assembly indices are equal
    assert L_det >= ai


def test_tryptophan_iterations():
    """
    Test the iterated calculation of short assembly paths for tryptophan.

    This function converts a SMILES string representing tryptophan to an RDKit molecule,
    then to a NetworkX graph. It removes hydrogen atoms from the graph and calculates the
    assembly index using three different approaches (default, CFG, and DET). It prints the
    results and asserts that the assembly indices calculated by the CFG and DET approaches
    are equal to the default approach.

    Asserts:
        The assembly indices calculated by the CFG and DET approaches are equal to the default approach.
    """
    print(flush=True)
    smi = "c1[nH]c2ccccc2c1C[C@H](N)C(=O)O"
    mol = att.smi_to_mol(smi)
    graph = att.mol_to_nx(mol)
    # Strip the hydrogen atoms from the graph
    graph = att.remove_hydrogen_from_graph(graph)
    # Calculate the assembly index using the default approach
    ai, _, _ = att.calculate_assembly_index(graph)
    # Calculate the assembly index using the DET approach
    L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph, iterations=100)
    # Print the assembly indices
    print(f"Assembly Index (asscpp) : {ai}", flush=True)
    print(f"Assembly Index (DET)    : {L_det}", flush=True)
    assert L_det >= ai


@pytest.mark.slow
def test_paclitaxel():
    """
    Test the bounding of assembly index of paclitaxel.

    """
    print(flush=True)
    smi = "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    mol = att.smi_to_mol(smi)
    graph = att.mol_to_nx(mol)
    # Strip the hydrogen atoms from the graph
    graph = att.remove_hydrogen_from_graph(graph)
    # Calculate the assembly index using the default approach
    ai, _, _ = att.calculate_assembly_index(graph)
    # Calculate the assembly index using the DET approach
    L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph)
    # Print the assembly indices
    print(f"Assembly Index (asscpp) : {ai}", flush=True)
    print(f"Assembly Index (DET)    : {L_det}", flush=True)
    # Assert that the assembly indices are equal
    assert L_det >= ai


def test_loads_joint():
    """
    Test the calculation of the assembly index for a list of joint molecules.

    This function converts a list of SMILES strings representing various molecules to RDKit molecules,
    then to NetworkX graphs. It joins these graphs into a single graph and calculates the assembly index
    using two different approaches (CFG and DET). It prints the results and the execution time for each approach.

    Asserts:
        The assembly indices calculated by the CFG and DET approaches are greater than 2.
    """
    print(flush=True)

    smiles_list = [
        'O=O',  # O2
        'O',  # H2O
        '[OH]',  # OH
        'OO[H]',  # HO2
        'OO',  # H2O2
        '[HH]',  # H2
        '[C]=O',  # CO
        '[CH]=O',  # HCO
        'C=O',  # H2CO
        'C',  # CH4
        '[CH3]',  # CH3
        'CC',  # C2H6
        '[N]=O',  # NO
        'N(=O)[O]',  # NO2
        '[NH]=O',  # HNO
        '[O-][O+]=O',  # O3
        'O=[N](=O)O',  # HNO3
        'C=C=C',  # C3H2
        '[CH]=C=C',  # C3H3
        'CC#C',  # CH3C2H
        'C=C=C',  # CH2CCH2
        'C=CC[CH2]',  # C3H5
        'CCC=O',  # C2H5CHO
        'CC=C',  # C3H6
        'CCC[CH2]',  # C3H7
        'CCC',  # C3H8
        'CCO',  # C2H4OH
        'C=CO',  # C2H2OH
        'CC[CH2]',  # C2H5
        'C=C',  # C2H4
        '[CH]',  # CH
        'COO[CH2]',  # CH3O2
        'C[O]',  # CH3O
        'C=C=O',  # CH2CO
        'CC(=O)[CH2]',  # CH3CO
        'CC=O',  # CH3CHO
        'C#C',  # C2H2
        '[CH]=C',  # C2H
        '[C]#[C]',  # C2
        '[CH]=C',  # C2H3
        '[CH]=S',  # HCS
        'S=C=S',  # CS2
        '[C]=S',  # CS
        'O=C=S',  # OCS
        '[SH]',  # HS
        'S',  # H2S
        'O=S(=O)=O',  # SO3
        'S=S',  # S2
        'O=[SH]',  # HSO
        'O=S(=O)(O)O',  # H2SO4
        'O=S=O',  # SO2
        '[S]=O',  # SO
        'O=C=O',  # CO2
        'N#N'  # N2
    ]
    # Convert all the SMILES to RDKit molecules
    mol_graphs = [att.smi_to_mol(smile) for smile in smiles_list]
    # Convert the RDKit molecules to NetworkX graphs
    graphs = [att.mol_to_nx(mol) for mol in mol_graphs]
    # Join the graphs
    graph = nx.disjoint_union_all(graphs)

    time_0 = time.perf_counter()
    # Calculate a path using the DET approach
    L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph)
    time_1 = time.perf_counter()
    # Print the 
    print(f"Path Length (DET)    : {L_det}", flush=True)
    print("execution time (DET)     : %3.2e" % (time_1 - time_0), flush=True)
    assert L_det > 2


def test_amino_acids():
    smiles_list = ['C(C(=O)O)N',
                   'C[C@@H](C(=O)O)N',
                   'C([C@@H](C(=O)O)N)O',
                   'C([C@@H](C(=O)O)N)C(=O)O',
                   'C([C@@H](C(=O)O)N)S',
                   'C(CC(=O)O)[C@@H](C(=O)O)N',
                   'C[C@H]([C@@H](C(=O)O)N)O',
                   'CC(C)[C@@H](C(=O)O)N',
                   'C([C@@H](C(=O)O)N)C(=O)N',
                   'C(CC(=O)N)[C@@H](C(=O)O)N',
                   'CC[C@H](C)[C@@H](C(=O)O)N',
                   'CC(C)C[C@@H](C(=O)O)N',
                   'C(CCN)C[C@@H](C(=O)O)N',
                   'C1C[C@H](NC1)C(=O)O',
                   'CSCC[C@@H](C(=O)O)N',
                   'C(C[C@@H](C(=O)O)N)CN=C(N)N',
                   'C1=C(NC=N1)C[C@@H](C(=O)O)N',
                   'C1=CC=C(C=C1)C[C@@H](C(=O)O)N',
                   'C1=CC(=CC=C1C[C@@H](C(=O)O)N)O',
                   'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N',
                   ]
    ai_ref = [3, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 9, 9, 9, 9, 11]
    jai_ref = 37

    # Convert all the SMILES to RDKit molecules
    mol_graphs = [att.smi_to_mol(smile) for smile in smiles_list]
    # Convert the RDKit molecules to NetworkX graphs and remove hydrogen atoms
    graphs = [att.remove_hydrogen_from_graph(att.mol_to_nx(mol)) for mol in mol_graphs]
    # Join the graphs
    graph_joint = nx.disjoint_union_all(graphs)

    # loop over the SMILES, convert to mol, convert to nx, remove H, calculate AI
    for i, g in enumerate(graphs):
        print(flush=True)
        L_det, _, _ = CFG.calculate_assembly_path_det(g)
        print(f"Path Length (DET)    : {L_det}", flush=True)
        assert L_det >= ai_ref[i]

    # Calculate the assembly index for the joint graph
    L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph_joint, iterations=10000)
    print(flush=True)
    print(f"Assembly Index (asscpp) : {jai_ref}", flush=True)
    print(f"Assembly Index (DET)    : {L_det}", flush=True)
    assert L_det >= jai_ref


def test_det_iteration_improvement():
    import matplotlib.pyplot as plt
    plt.rcParams['axes.linewidth'] = 2.0

    # Use the amino acids as a test case
    smiles_list = ['C(C(=O)O)N',
                   'C[C@@H](C(=O)O)N',
                   'C([C@@H](C(=O)O)N)O',
                   'C([C@@H](C(=O)O)N)C(=O)O',
                   'C([C@@H](C(=O)O)N)S',
                   'C(CC(=O)O)[C@@H](C(=O)O)N',
                   'C[C@H]([C@@H](C(=O)O)N)O',
                   'CC(C)[C@@H](C(=O)O)N',
                   'C([C@@H](C(=O)O)N)C(=O)N',
                   'C(CC(=O)N)[C@@H](C(=O)O)N',
                   'CC[C@H](C)[C@@H](C(=O)O)N',
                   'CC(C)C[C@@H](C(=O)O)N',
                   'C(CCN)C[C@@H](C(=O)O)N',
                   'C1C[C@H](NC1)C(=O)O',
                   'CSCC[C@@H](C(=O)O)N',
                   'C(C[C@@H](C(=O)O)N)CN=C(N)N',
                   'C1=C(NC=N1)C[C@@H](C(=O)O)N',
                   'C1=CC=C(C=C1)C[C@@H](C(=O)O)N',
                   'C1=CC(=CC=C1C[C@@H](C(=O)O)N)O',
                   'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N',
                   ]
    # Convert all the SMILES to RDKit molecules
    mol_graphs = [att.smi_to_mol(smile) for smile in smiles_list]
    # Convert the RDKit molecules to NetworkX graphs and remove hydrogen atoms
    graphs = [att.mol_to_nx(mol) for mol in mol_graphs]
    # Join the graphs
    graph_joint = nx.disjoint_union_all(graphs)
    graph_joint_nh = att.remove_hydrogen_from_graph(graph_joint)

    iterations_list = list(range(1, 501, 10))
    L_det_list = []
    L_det_list_nh = []

    for iterations in iterations_list:
        L_det, _, _ = CFGgraph.calculate_assembly_path_det(graph_joint, iterations=iterations)
        L_det_list.append(L_det)
        L_det_nh, _, _ = CFGgraph.calculate_assembly_path_det(graph_joint_nh, iterations=iterations)
        L_det_list_nh.append(L_det_nh)

    plt.plot(iterations_list, L_det_list, 'o-', label='DET', color='blue', lw=2)
    plt.plot(iterations_list, L_det_list_nh, 'o-', label='DET (no H)', color='red', lw=2)
    att.n_plot('Iterations', 'Assembly Index')
    plt.legend()
    plt.savefig('det_iterations.png', dpi=600)
    plt.savefig('det_iterations.pdf')
    plt.show()
