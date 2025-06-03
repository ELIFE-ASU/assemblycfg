import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from typing import List


def safe_standardize_mol(mol: Chem.Mol, add_hydrogens: bool = True) -> Chem.Mol:
    """
    Standardise the given RDKit molecule with additional safety checks.

    Args:
        mol (rdkit.Chem.Mol): The input RDKit molecule to be standardised.
        add_hydrogens (bool, optional): Whether to add hydrogens to the molecule. Default is True.

    Returns:
        rdkit.Chem.Mol: The standardised RDKit molecule.
    """
    # Update the molecule's property cache without strict checking
    mol.UpdatePropertyCache(strict=False)
    # Set conjugation and hybridisation states
    Chem.SetConjugation(mol)
    Chem.SetHybridization(mol)
    # Normalise the molecule, excluding clean-up and property sanitisation
    Chem.SanitizeMol(mol, sanitizeOps=(Chem.SANITIZE_ALL ^ Chem.SANITIZE_CLEANUP ^ Chem.SANITIZE_PROPERTIES),
                     catchErrors=False)
    # Normalise the molecule in place using RDKit's MolStandardize
    rdMolStandardize.NormalizeInPlace(mol)
    # Kekulise the molecule (convert aromatic bonds to alternating single and double bonds)
    Chem.Kekulize(mol)
    if add_hydrogens:
        # Add hydrogens
        mol = Chem.AddHs(mol)
    return mol


def mol_to_nx(mol: Chem.Mol, add_hydrogens: bool = True) -> nx.Graph:
    """
    Convert an RDKit molecule to a NetworkX graph.

    Args:
        mol (Chem.Mol): The RDKit molecule to convert.
        add_hydrogens (bool, optional): Whether to keep hydrogen atoms in the graph. Default is True.

    Returns:
        nx.Graph: The resulting NetworkX graph where nodes represent atoms and edges represent bonds.
    """

    mol = safe_standardize_mol(mol, add_hydrogens=add_hydrogens)

    graph = nx.Graph()

    for atom in mol.GetAtoms():
        graph.add_node(atom.GetIdx(), color=atom.GetSymbol())

    for bond in mol.GetBonds():
        bond_type = bond_order_rdkit_to_int(bond.GetBondType())
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), color=bond_type)

    if not add_hydrogens:
        graph = remove_hydrogen_from_graph(graph)

    return graph


def remove_hydrogen_from_graph(graph: nx.Graph) -> nx.Graph:
    """
    Remove all hydrogen atoms from a NetworkX graph.

    Args:
        graph (nx.Graph): The input NetworkX graph where nodes represent atoms.

    Returns:
        nx.Graph: The modified graph with all hydrogen atoms removed.
    """
    nodes = list(graph.nodes())
    for node in nodes:
        if graph.nodes[node]["color"] == "H":
            graph.remove_node(node)
    return graph


def bond_order_rdkit_to_int(bond_type: Chem.BondType) -> int:
    """
    Convert RDKit's BondType to a bond order int.

    Args:
        bond_type (Chem.BondType): The RDKit BondType to convert.

    Returns:
        int: The corresponding bond order int.
    """
    converter = {
        Chem.rdchem.BondType.SINGLE: 1,
        Chem.rdchem.BondType.DOUBLE: 2,
        Chem.rdchem.BondType.TRIPLE: 3,
        Chem.rdchem.BondType.QUADRUPLE: 4,
        Chem.rdchem.BondType.QUINTUPLE: 5,
        Chem.rdchem.BondType.IONIC: 6
    }
    return converter.get(bond_type, 1)


def get_disconnected_subgraphs(graph: nx.Graph) -> List[nx.Graph]:
    """
    Return subgraphs of connected components without copying if not necessary.

    Args:
        graph (nx.Graph): The input graph.

    Returns:
        List[nx.Graph]: A list of subgraphs, each representing a connected component.
    """
    return [graph.subgraph(c) for c in nx.connected_components(graph)]


def molfile_to_mol(mol: str, add_hydrogens: bool = True, safe_sanitise: bool = False) -> Chem.Mol:
    """
    Convert a Molfile to a standardized RDKit molecule.

    Args:
        mol (str): The path to the Molfile representing the molecule.
        add_hydrogens (bool, optional): Whether to add hydrogens to the molecule. Default is True.
        safe_sanitise (bool, optional): Whether to use safe sanitisation. Default is False.

    Returns:
        Chem.Mol: The standardized RDKit molecule.
    """
    # Convert the Molfile to an RDKit molecule
    mol = Chem.MolFromMolFile(mol)
    # Standardise the molecule
    if safe_sanitise:
        return safe_standardize_mol(mol, add_hydrogens=add_hydrogens)
    else:
        return standardize_mol(mol, add_hydrogens=add_hydrogens)


def standardize_mol(mol: Chem.Mol, add_hydrogens: bool = True) -> Chem.Mol:
    """
    Standardize the given RDKit molecule.

    Args:
        mol (Chem.Mol): The input RDKit molecule to be standardised.
        add_hydrogens (bool, optional): Whether to add hydrogens to the molecule. Default is True.

    Returns:
        Chem.Mol: The standardized RDKit molecule.
    """
    # Sanitise the molecule
    Chem.SanitizeMol(mol, catchErrors=False)
    # Normalise the molecule in place using RDKit's MolStandardize
    rdMolStandardize.NormalizeInPlace(mol)
    # Update the molecule's property cache without strict checking
    mol.UpdatePropertyCache(strict=True)
    # Kekulise the molecule (convert aromatic bonds to alternating single and double bonds)
    Chem.Kekulize(mol)
    if add_hydrogens:
        # Add hydrogens
        mol = Chem.AddHs(mol)
    # Return the molecule
    return mol


def nx_to_dict(graph):
    """
    Convert a NetworkX graph to a dictionary with node colors, neighbors, and edge colors.

    Args:
        graph (nx.Graph): The input NetworkX graph.

    Returns:
        dict: A dictionary representation of the graph with node colors, neighbors, and edge colors.
    """
    graph_dict = {}
    for node_id, node_data in graph.nodes(data=True):
        neighbors = list(graph.neighbors(node_id))
        edge_colors = [graph[node_id][neighbor].get('color', 1) for neighbor in neighbors]
        graph_dict[node_id] = {
            'vertex_color': node_data.get('color', 'C'),  # Default color is 'C' if not specified
            'neighbors': neighbors,
            'edge_colors': edge_colors
        }
    return graph_dict


def dict_to_nx(graph_dict):
    """
    Convert a dictionary representation of a graph back to a NetworkX graph.

    Args:
        graph_dict (dict): A dictionary representation of the graph with node colors, neighbors, and edge colors.

    Returns:
        nx.Graph: The reconstructed NetworkX graph.
    """
    graph = nx.Graph()

    # Add nodes with their attributes
    for node_id, node_data in graph_dict.items():
        graph.add_node(node_id, color=node_data['vertex_color'])

    # Add edges with their attributes
    for node_id, node_data in graph_dict.items():
        for neighbor, edge_color in zip(node_data['neighbors'], node_data['edge_colors']):
            if not graph.has_edge(node_id, neighbor):
                graph.add_edge(node_id, neighbor, color=edge_color)

    return graph


def mol2graph(mol_file):
    """
    Convert a molecular file to a graph dictionary representation.

    Args:
        mol_file (str): The path to the molecular file.

    Returns:
        dict: A dictionary representation of the molecular graph with node colors, neighbors, and edge colors.
    """
    mol = molfile_to_mol(mol_file)  # Convert the molecular file to a molecule object
    mol_graph = mol_to_nx(mol)  # Convert the molecule object to a NetworkX graph
    mol_graph = remove_hydrogen_from_graph(mol_graph)  # Remove hydrogen atoms from the graph
    return nx_to_dict(mol_graph)  # Convert the NetworkX graph to a dictionary representation


def print_graph_dict(graph_dict):
    """
    Print the contents of a graph dictionary.

    Args:
        graph_dict (dict): A dictionary representation of a graph with node colors, neighbors, and edge colors.
    """
    for i in graph_dict:
        print(f"{i}: {graph_dict[i]}", flush=True)
    print(flush=True)


def print_virtual_objects(virtual_objects):
    """
    Print the contents of a list of virtual objects, where each virtual object is a NetworkX graph.

    Args:
        virtual_objects (list of nx.Graph): A list of NetworkX graph objects representing virtual objects.
    """
    for i, item in enumerate(virtual_objects):
        print(f"{i}: ", flush=True)
        print_graph_dict(nx_to_dict(item))


def nx_to_dict_v2(graph):
    """
    Convert a NetworkX graph to a dictionary with vertex colors and edges with bond types.

    Args:
        graph (nx.Graph): The input NetworkX graph.

    Returns:
        tuple: A tuple containing:
            - vmap (dict): A dictionary mapping nodes to their colors.
            - edges (list): A list of tuples representing the edges with bond types, where each tuple is (node1, node2, bond_type).
    """
    # Create vertex map
    vmap = {node: data['color'] for node, data in graph.nodes(data=True)}

    # Create edges list with bond types
    edges = []
    for u, v, data in graph.edges(data=True):
        bond_type = f"B{data['color']}"
        edges.append((u, v, bond_type))

    return vmap, edges


def dict_to_nx_v2(graph_dict):
    """
    Convert a dictionary representation of a graph back to a NetworkX graph.

    Args:
        graph_dict (tuple): A tuple containing:
            - vmap (dict): A dictionary mapping nodes to their colors.
            - edges (list): A list of tuples representing the edges with bond types, where each tuple is (node1, node2, bond_type).

    Returns:
        nx.Graph: The reconstructed NetworkX graph.
    """
    vmap, edges = graph_dict
    graph = nx.Graph()

    # Add nodes with their attributes
    for node, color in vmap.items():
        graph.add_node(node, color=color)

    # Add edges with their attributes
    for u, v, bond_type in edges:
        bond_color = int(bond_type[1:])  # Extract the integer part of the bond type
        graph.add_edge(u, v, color=bond_color)

    return graph
