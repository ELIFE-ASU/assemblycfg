import time

import networkx as nx

import CFG


def test_basic():
    """
    Test the calculation of a short assembly path using the DET method. 
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
    graph = CFG.dict_to_nx(d)

    # Calculate a short path using DET
    L_det, _, _ = CFG.calculate_assembly_path_det(graph)

    print(f"Path Length (DET)    : {L_det}", flush=True)

    # Assert that DET finds a path no shorter than the assembly index
    assert L_det >= 3


def test_basic_joint():
    """
    Test the calculation of a short assembly path for a joint system of molecules.
    """
    print(flush=True)
    molecules = "O=C=O.C.O.N"
    molecules = molecules.split(".")
    # Convert all the SMILES to networkX graphs
    graphs = [CFG.smi_to_nx(smile) for smile in molecules]
    # Join the graphs into a single graph
    graph = nx.disjoint_union_all(graphs)

    # Calculate the assembly index using the DET approach
    L_det, _, _ = CFG.calculate_assembly_path_det(graph)
    print(f"Path Length (DET)    : {L_det}", flush=True)

    # Assert that DET finds a path no shorter than the assembly index
    assert L_det >= 6


def test_det_cyclobutane():
    """
    Test the calculation of the assembly index for cyclobutane.
    """
    print(flush=True)
    smi = "C1CCC1"
    graph = CFG.smi_to_nx(smi)

    L_det, _, _ = CFG.calculate_assembly_path_det(graph)
    print(f"Path Length (DET)    : {L_det}", flush=True)

    print("Without hydrogen...", flush=True)
    graph = CFG.smi_to_nx(smi, add_hydrogens=False)

    # Calculate the assembly index using the DET approach
    L_det_nh, _, _ = CFG.calculate_assembly_path_det(graph)
    print(f"Path Length (DET)    : {L_det_nh}", flush=True)

    # Assert that the paths are no shorter than the assembly indices
    assert L_det >= 4
    assert L_det_nh >= 2


def test_long_boi():
    """
    Test the calculation of a short assembly path for a long alkane chain.
    """
    print(flush=True)
    smi = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    graph = CFG.smi_to_nx(smi, add_hydrogens=False)

    L_det, _, _ = CFG.calculate_assembly_path_det(graph)

    print(f"Path Length (DET)    : {L_det}", flush=True)

    # Assert that the path is no shorter than the assembly index
    assert L_det >= 7


def test_benzene():
    """
    Test DET on benzene.
    """
    print(flush=True)
    smi = "c1ccccc1"
    graph = CFG.smi_to_nx(smi)
    L_det, _, _ = CFG.calculate_assembly_path_det(graph)
    print(f"Path Length (DET)    : {L_det}", flush=True)
    # Assert that the path is no shorter than the assembly index
    assert L_det >= 5


def test_tryptophan():
    """
    Test det on tryptophan.
    """
    print(flush=True)
    smi = "c1[nH]c2ccccc2c1C[C@H](N)C(=O)O"
    graph = CFG.smi_to_nx(smi, add_hydrogens=False)
    L_det, _, _ = CFG.calculate_assembly_path_det(graph)

    print(f"Path Length (DET)    : {L_det}", flush=True)
    # Assert that the path is no shorter than the assembly index
    assert L_det >= 11


def test_tryptophan_iterations():
    """
    Test the iterated calculation of short assembly paths for tryptophan.
    """
    print(flush=True)
    smi = "c1[nH]c2ccccc2c1C[C@H](N)C(=O)O"
    graph = CFG.smi_to_nx(smi, add_hydrogens=False)
    L_det, _, _ = CFG.calculate_assembly_path_det(graph, iterations=100)

    print(f"Path Length (DET)    : {L_det}", flush=True)
    assert L_det >= 11  # Assert that the path is no shorter than the assembly index


def test_loads_joint():
    """
    Test det a list of molecules.
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
    # Convert all the SMILES to NetworkX graphs
    graphs = [CFG.smi_to_nx(smi) for smi in smiles_list]
    # Join the graphs
    graph = nx.disjoint_union_all(graphs)

    time_0 = time.perf_counter()
    # Calculate a path using the DET approach
    L_det, _, _ = CFG.calculate_assembly_path_det(graph)
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

    # Convert all the SMILES to NetworkX graphs
    graphs = [CFG.smi_to_nx(smile, add_hydrogens=False) for smile in smiles_list]
    # Join the graphs
    graph_joint = nx.disjoint_union_all(graphs)

    # loop over the graphs and calculate a short assembly path using DET
    for i, g in enumerate(graphs):
        print(flush=True)
        L_det, _, _ = CFG.calculate_assembly_path_det(g)
        print(f"Path Length (DET)    : {L_det}", flush=True)
        assert L_det >= ai_ref[i]  # Assert that the path is no shorter than the assembly index

    # Calculate the assembly index for the joint graph
    L_det, _, _ = CFG.calculate_assembly_path_det(graph_joint, iterations=1000)
    print(flush=True)
    print(f"Assembly Index (asscpp) : {jai_ref}", flush=True)
    print(f"Path Length (DET)    : {L_det}", flush=True)
    assert L_det >= jai_ref  # Assert that the path is no shorter than the assembly index
