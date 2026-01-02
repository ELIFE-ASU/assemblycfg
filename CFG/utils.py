import warnings
from typing import List, Dict, Any, Sequence, Optional

import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def safe_standardize_mol(mol: Chem.Mol, add_hydrogens: bool = True) -> Chem.Mol:
    """
    Standardize an RDKit molecule with relaxed safety checks and optional hydrogen addition.

    Performs a sequence of RDKit normalisation and sanitisation steps intended to
    be more tolerant of imperfect input molecules: property-cache update without
    strict checking, explicit conjugation/hybridisation setting, selective
    sanitisation that omits cleanup/property operations, in-place normalisation
    via MolStandardize, kekulisation, and optional hydrogen addition.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Input RDKit molecule to standardize. The function operates in-place on
        the provided molecule object for most operations but may return a new
        molecule instance when hydrogens are added.
    add_hydrogens : bool, optional
        If True, explicit hydrogens are added to the molecule using
        ``Chem.AddHs`` after sanitisation and kekulisation. Default is True.

    Returns
    -------
    rdkit.Chem.Mol
        The standardized RDKit molecule. If ``add_hydrogens`` is True a new
        molecule instance that includes explicit hydrogens may be returned;
        otherwise the original (mutated) molecule is returned.
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


def smi_to_mol(smi: str, add_hydrogens: bool = True, sanitize: bool = True) -> Optional[Chem.Mol]:
    """
    Convert a SMILES string to an RDKit molecule object.

    This function takes a SMILES (Simplified Molecular Input Line Entry System) string
    and converts it into an RDKit molecule object. Optionally, the molecule can be
    sanitized and explicit hydrogens can be added.

    Parameters:
    -----------
    smi : str
        A SMILES string representing the molecular structure.
    add_hydrogens : bool, optional
        If True, adds explicit hydrogens to the molecule during sanitization. Defaults to True.
    sanitize : bool, optional
        If True, sanitizes the molecule after conversion. Defaults to True.

    Returns:
    --------
    Chem.Mol
        An RDKit molecule object representing the input SMILES string.

    Raises:
    -------
    Warning
        If the SMILES string contains disconnected molecules (indicated by a '.' character).

    Notes:
    ------
    - The function uses RDKit's `MolFromSmiles` to create the molecule object.
    - Sanitization ensures the molecule is chemically valid and standardized.
    - Disconnected molecules in the SMILES string are flagged with a warning.
    """
    if '.' in smi:
        warnings.warn("Disconnected molecules detected in SMILES string. Ensure proper handling of these molecules.")
    # Convert the SMILES string to an RDKit molecule object
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    # Sanitise the molecule
    if sanitize:
        return safe_standardize_mol(mol, add_hydrogens=add_hydrogens)
    else:
        return mol


def smi_to_nx(smi: str, add_hydrogens: bool = True) -> nx.Graph:
    """
    Convert a SMILES string to a NetworkX molecular graph.

    Reads a SMILES string, converts it to an RDKit ``Mol`` using
    ``smi_to_mol`` (which performs optional sanitisation and hydrogen
    addition), then converts the resulting molecule into a NetworkX graph
    via ``mol_to_nx``.

    Parameters
    ----------
    smi : str
        SMILES string describing the molecular structure.
    add_hydrogens : bool, optional
        If True (default), explicit hydrogens produced during molecule
        standardisation are preserved in the resulting graph. If False,
        hydrogen nodes will be removed from the graph.

    Returns
    -------
    networkx.Graph
        Undirected NetworkX graph where each node id is the RDKit atom
        index and each node has a ``'color'`` attribute (atomic symbol).
        Edges represent bonds and carry an integer ``'color'`` attribute
        encoding bond order.

    Raises
    ------
    TypeError
        If ``smi`` is not a ``str``.
    ValueError
        If RDKit fails to parse the SMILES into a valid molecule (for
        example when ``smi_to_mol`` returns ``None``).
    """
    mol = smi_to_mol(smi, add_hydrogens=add_hydrogens)
    return mol_to_nx(mol, add_hydrogens=add_hydrogens)


def mol_to_nx(mol: Chem.Mol, add_hydrogens: bool = True) -> nx.Graph:
    """
    Convert an RDKit molecule to a NetworkX graph.

    Constructs an undirected NetworkX graph from an RDKit ``Mol`` by first
    standardizing the molecule (via ``safe_standardize_mol``), creating a node
    for each atom (node id = RDKit atom index) with a ``'color'`` node attribute
    set to the atom symbol, and adding edges for each bond with an integer
    ``'color'`` attribute representing bond order.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule to convert. The function will call ``safe_standardize_mol``
        on this object; callers should pass a valid RDKit ``Mol`` instance.
    add_hydrogens : bool, optional
        If True (default), hydrogens produced by standardisation are kept in the
        resulting graph. If False, hydrogen nodes are removed after graph
        construction via ``remove_hydrogen_from_graph``.

    Returns
    -------
    networkx.Graph
        Undirected graph where node identifiers are RDKit atom indices (integers)
        and each node has a ``'color'`` attribute (atom symbol). Edges correspond
        to bonds and carry a ``'color'`` attribute containing the bond order as
        an integer (converted with ``bond_order_rdkit_to_int``).

    Raises
    ------
    TypeError
        If ``mol`` is not an RDKit ``Mol`` instance.
    ValueError
        If molecule parsing or standardisation fails (for example when RDKit
        returns ``None`` for the input).
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
    Remove all hydrogen atoms from a NetworkX molecular graph.

    Removes nodes whose ``'color'`` node attribute equals ``'H'``. The operation
    is performed by iterating over a static list of node identifiers to allow
    in-place removal without iterator invalidation.

    Parameters
    ----------
    graph : networkx.Graph
        Input NetworkX graph where nodes represent atoms and each node has a
        ``'color'`` attribute containing the atomic symbol (e.g. ``'C'``, ``'H'``).

    Returns
    -------
    networkx.Graph
        The same graph instance with all hydrogen nodes removed. The function
        mutates the provided graph and also returns it for convenience.

    Raises
    ------
    TypeError
        If ``graph`` is not a ``networkx.Graph`` instance.
    KeyError
        If a node in ``graph`` does not have a ``'color'`` attribute.
    """
    nodes = list(graph.nodes())
    for node in nodes:
        if graph.nodes[node]["color"] == "H":
            graph.remove_node(node)
    return graph


def bond_order_rdkit_to_int(bond_type: Chem.BondType) -> int:
    """
    Convert an RDKit bond type to an integer bond order.

    Maps common RDKit ``BondType`` values to integer bond orders used in
    graph-based representations. Unknown or unhandled bond types return a
    conservative default of ``1``.

    Parameters
    ----------
    bond_type : rdkit.Chem.rdchem.BondType or Chem.BondType
        RDKit bond type object (for example as returned by ``bond.GetBondType()``).

    Returns
    -------
    int
        Integer bond order corresponding to ``bond_type``. Typical mappings:
        - ``BondType.SINGLE`` -> 1
        - ``BondType.DOUBLE`` -> 2
        - ``BondType.TRIPLE`` -> 3
        - ``BondType.QUADRUPLE`` -> 4
        - ``BondType.QUINTUPLE`` -> 5
        - ``BondType.IONIC`` -> 6
        If ``bond_type`` is not recognized the function returns ``1`` as a default.

    Raises
    ------
    None
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
    Return subgraphs for each connected component of a NetworkX graph.

    Constructs a list of subgraph views corresponding to every connected component
    in the input graph. The returned subgraphs are views into the original graph
    when possible (they are not deep copies); call ``.copy()`` on each subgraph if
    an independent copy is required.

    Parameters
    ----------
    graph : networkx.Graph
        Input undirected graph whose connected components will be extracted.

    Returns
    -------
    List[nx.Graph]
        A list of subgraph objects (one per connected component). Each subgraph
        contains the nodes and edges that belong to that component. Modifications
        to the original graph may affect these subgraph views.

    Raises
    ------
    TypeError
        If ``graph`` is not a ``networkx.Graph`` instance.
    """
    return [graph.subgraph(c) for c in nx.connected_components(graph)]


def molfile_to_mol(mol: str, add_hydrogens: bool = True, safe_sanitise: bool = False) -> Optional[Chem.Mol]:
    """
    Convert a Molfile (file path) to a standardized RDKit molecule.

    Reads a Molfile using RDKit and returns a standardized RDKit ``Chem.Mol``.
    Standardisation is performed by either ``safe_standardize_mol`` (when
    ``safe_sanitise`` is True) or ``standardize_mol`` (when False).

    Parameters
    ----------
    mol : str
        Path to the Molfile to read. The function calls ``Chem.MolFromMolFile`` on
        this path; passing raw Molfile content is not supported by this function.
    add_hydrogens : bool, optional
        If True, explicit hydrogens will be added to the returned molecule.
        Default is True.
    safe_sanitise : bool, optional
        If True, perform relaxed sanitisation via ``safe_standardize_mol``.
        If False, use ``standardize_mol``. Default is False.

    Returns
    -------
    rdkit.Chem.Mol or None
        A standardized RDKit ``Mol`` instance on success. If RDKit fails to parse
        the input file the function may return ``None``; callers should check the
        return value before use.

    Raises
    ------
    TypeError
        If ``mol`` is not a string.
    FileNotFoundError
        If the provided file path does not exist.
    ValueError
        If RDKit fails to produce a molecule from the given Molfile.
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

    Performs strict sanitisation and canonicalisation steps on an RDKit
    ``Chem.Mol``: sanitisation, in-place normalisation using
    ``rdMolStandardize.NormalizeInPlace``, property-cache update with strict
    checking, kekulisation, and optional explicit hydrogen addition.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Input RDKit molecule to be standardised. The function mutates the
        provided molecule for most operations and may return a new instance if
        hydrogens are added.
    add_hydrogens : bool, optional
        If True (default), explicit hydrogens are added to the returned
        molecule via ``Chem.AddHs``.

    Returns
    -------
    rdkit.Chem.Mol
        The standardized RDKit molecule. When ``add_hydrogens`` is True a new
        molecule instance containing explicit hydrogens may be returned;
        otherwise the original (mutated) molecule is returned.

    Raises
    ------
    TypeError
        If ``mol`` is not an instance of ``rdkit.Chem.Mol``.
    ValueError
        If RDKit sanitisation or parsing fails (for example when RDKit returns
        ``None`` or raises during sanitisation/kekulisation).
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


def nx_to_dict(graph: nx.Graph) -> Dict[Any, Dict[str, Any]]:
    """
    Convert a NetworkX graph to a dictionary with node colors, neighbors, and edge colors.

    Parameters
    ----------
    graph : networkx.Graph
        Input NetworkX graph. Nodes are expected to have a ``'color'`` attribute;
        edges may have a ``'color'`` attribute representing bond/order. If a node
        or edge attribute is missing sensible defaults are used (see Notes).

    Returns
    -------
    dict
        Mapping from node id to a dictionary with the following keys:
        - ``'vertex_color'`` : str
            Node color (e.g. atomic symbol). Defaults to ``'C'`` when the node
            does not define a ``'color'`` attribute.
        - ``'neighbors'`` : list[int]
            List of neighboring node identifiers.
        - ``'edge_colors'`` : list[int]
            List of edge color values (integers) aligned with ``'neighbors'``;
            each entry corresponds to the edge color for the neighbor at the same index.

    Raises
    ------
    TypeError
        If ``graph`` is not an instance of ``networkx.Graph``.
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


def dict_to_nx(graph_dict: Dict[Any, Dict[str, Any]]) -> nx.Graph:
    """
    Convert a dictionary representation of a graph back to a NetworkX graph.

    The input dictionary is expected to map node identifiers to dictionaries
    containing a ``'vertex_color'`` (string), a ``'neighbors'`` list of node ids,
    and an ``'edge_colors'`` list of integers aligned with ``'neighbors'``.

    Parameters
    ----------
    graph_dict : dict
        Mapping from node id to a dictionary with keys:
        - ``'vertex_color'`` (str): Node color (e.g. atomic symbol).
        - ``'neighbors'`` (list[int]): Neighbor node identifiers.
        - ``'edge_colors'`` (list[int]): Edge color values aligned with ``'neighbors'``.

    Returns
    -------
    networkx.Graph
        A NetworkX undirected graph reconstructed from ``graph_dict``. Each node
        is added with a ``'color'`` attribute set from ``'vertex_color'`` and each
        undirected edge is added once with a ``'color'`` attribute taken from the
        corresponding ``'edge_colors'`` entry.

    Raises
    ------
    TypeError
        If ``graph_dict`` is not a ``dict``.
    KeyError
        If a node entry does not contain the required keys
        (``'vertex_color'``, ``'neighbors'``, ``'edge_colors'``).
    ValueError
        If the lengths of ``'neighbors'`` and ``'edge_colors'`` for any node do
        not match.
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


def mol2graph(mol_file: str) -> Dict[Any, Dict[str, Any]]:
    """
    Convert a molecular file to a graph dictionary representation.

    Reads a molecular file using :func:`molfile_to_mol`, converts the resulting RDKit
    molecule to a NetworkX graph via :func:`mol_to_nx`, removes hydrogen atoms, and
    returns a standalone dictionary representation produced by :func:`nx_to_dict`.

    Parameters
    ----------
    mol_file : str
        Path to the molecular file to read (for example, `molecule.mol`).

    Returns
    -------
    dict
        Dictionary mapping node identifiers to dictionaries with keys:
        - ``'vertex_color'`` (str): atomic symbol for the node,
        - ``'neighbors'`` (list[int]): neighbor node identifiers,
        - ``'edge_colors'`` (list[int]): bond orders aligned with ``'neighbors'``.

    Raises
    ------
    TypeError
        If ``mol_file`` is not a string.
    FileNotFoundError
        If the provided file path does not exist.
    ValueError
        If RDKit fails to parse the file into a valid molecule.
    True
    """
    mol = molfile_to_mol(mol_file)  # Convert the molecular file to a molecule object
    mol_graph = mol_to_nx(mol)  # Convert the molecule object to a NetworkX graph
    mol_graph = remove_hydrogen_from_graph(mol_graph)  # Remove hydrogen atoms from the graph
    return nx_to_dict(mol_graph)  # Convert the NetworkX graph to a dictionary representation


def print_graph_dict(graph_dict: Dict[Any, Dict[str, Any]]) -> None:
    """
    Print the contents of a graph dictionary.

    Parameters
    ----------
    graph_dict : dict
        Dictionary mapping node identifiers to dictionaries with keys
        ``'vertex_color'``, ``'neighbors'``, and ``'edge_colors'`` (as produced by
        ``nx_to_dict``). Each node entry is printed on its own line to standard
        output using Python's ``print`` and flushed immediately.

    Returns
    -------
    None
        This function prints to stdout and returns ``None``.

    Raises
    ------
    TypeError
        If ``graph_dict`` is not a ``dict``.
    KeyError
        If a node entry does not contain the expected keys.
    """
    for i in graph_dict:
        print(f"{i}: {graph_dict[i]}", flush=True)
    print(flush=True)


def print_virtual_objects(virtual_objects: Sequence[nx.Graph]) -> None:
    """
    Print a sequence of virtual objects (NetworkX graphs) to standard output.

    Iterates over the provided sequence of NetworkX graph objects and prints each
    graph's dictionary representation using ``print_graph_dict``. Each graph is
    prefixed by its zero-based index on its own line.

    Parameters
    ----------
    virtual_objects : Sequence[nx.Graph]
        Iterable of NetworkX graphs representing virtual objects. Each item is
        converted with ``nx_to_dict`` before printing.

    Returns
    -------
    None
        This function prints to stdout for human-readable inspection and returns
        ``None``.

    Raises
    ------
    TypeError
        If ``virtual_objects`` is not iterable or contains items that are not
        instances of ``networkx.Graph``.
    """
    for i, item in enumerate(virtual_objects):
        print(f"{i}: ", flush=True)
        print_graph_dict(nx_to_dict(item))
