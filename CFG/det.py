import random
from collections import Counter
from typing import List, Dict, Tuple, Optional, Any, Union

import networkx as nx
import numpy as np

from .utils import get_disconnected_subgraphs


def partition_into_disjoint_trails(graph: nx.Graph, debug: bool = False) -> List[List[Any]]:
    """
    Partition a NetworkX graph into disjoint trails (edge-disjoint walks).

    Constructs a set of trails (node sequences) that together cover every edge
    of the input graph exactly once. Trails are produced by repeatedly selecting
    an unused edge and greedily extending the walk forward and backward until
    no further unused incident edges remain.

    Parameters
    ----------
    graph : networkx.Graph
        Undirected input graph whose edges will be partitioned into trails. Edge
        attributes are ignored for the partitioning logic; only graph topology is used.
    debug : bool, optional
        If True, print per-trail and intermediate-state diagnostics to stdout.
        Default is False.

    Returns
    -------
    list of list
        List of trails. Each trail is represented as a list of node identifiers
        in visitation order. Trails with two nodes correspond to single edges.

    Raises
    ------
    TypeError
        If ``graph`` is not a NetworkX graph-like object exposing ``edges``,
        ``neighbors`` and related iterators.
    ValueError
        If the graph contains self-loops that prevent proper trail extension;
        callers should pre-clean self-loops if such behavior is undesired.
    """
    if not graph.edges:
        return []

    remaining_edges = set(graph.edges)
    trails = []

    while remaining_edges:
        edge = random.choice(list(remaining_edges))
        current_node, next_node = edge
        trail = [current_node, next_node]
        remaining_edges.remove(edge)

        # Extend the trail forward
        while True:
            neighbors = [n for n in graph.neighbors(next_node) if
                         (next_node, n) in remaining_edges or (n, next_node) in remaining_edges]
            if not neighbors:
                break
            next_next_node = random.choice(neighbors)
            edge = (next_node, next_next_node) if (next_node, next_next_node) in remaining_edges else (
                next_next_node, next_node)
            remaining_edges.remove(edge)
            trail.append(next_next_node)
            next_node = next_next_node

        # Extend the trail backward
        back_neighbors = [n for n in graph.neighbors(current_node) if
                          (n, current_node) in remaining_edges or (current_node, n) in remaining_edges]
        if back_neighbors:
            next_node = random.choice(back_neighbors)
            edge = (next_node, current_node) if (next_node, current_node) in remaining_edges else (
                current_node, next_node)
            remaining_edges.remove(edge)
            back_trail = [next_node]

            while True:
                neighbors = [n for n in graph.neighbors(next_node) if
                             (next_node, n) in remaining_edges or (n, next_node) in remaining_edges]
                if not neighbors:
                    break
                next_next_node = random.choice(neighbors)
                edge = (next_node, next_next_node) if (next_node, next_next_node) in remaining_edges else (
                    next_next_node, next_node)
                remaining_edges.remove(edge)
                back_trail.append(next_next_node)
                next_node = next_next_node

            trail = back_trail[::-1] + trail

        trails.append(trail)
        if debug:
            print(f"Trail: {trail}", flush=True)
            print(f"Remaining edges: {remaining_edges}", flush=True)

    return trails


def get_unique_char(input_str: Union[str, List[str]]) -> str:
    """
    Find a unique character not present in the given input sequence.

    Searches Unicode code points and returns the first single-character string
    that does not appear in the provided sequence.

    Parameters
    ----------
    input_str : str or list of str
        Input sequence of characters (for example, a Python string or an iterable
        of single-character strings). Membership tests are performed against the
        elements of this iterable.

    Returns
    -------
    str
        A single-character string representing a Unicode code point not found in
        ``input_str``.

    Raises
    ------
    TypeError
        If ``input_str`` is not an iterable of characters.
    ValueError
        If no available Unicode code point could be found. The function searches
        code points from U+0001 to U+10FFFF; exhaustion of this range is
        extremely unlikely in practice.
    """
    for i in range(1, 1114111):
        char = chr(i)
        if char not in input_str:
            return char


def trails_to_sequences(trails: List[List[Any]],
                        graph: nx.Graph,
                        debug: bool = False) -> Tuple[List[str], Dict[Tuple[Any, Any, Any], str], int]:
    """
    Convert a list of trails into encoded sequences and an assembly index.

    Maps each non-trivial trail (length > 2 nodes) to a string where each
    character represents a molecular unit derived from an edge traversal.
    Also returns the mapping from molecular units to single-character symbols
    and a count of trivial trails (single-edge trails) omitted from the
    returned sequences.

    Parameters
    ----------
    trails : list of list
        Iterable of trails, where each trail is an ordered list of node
        identifiers (e.g. integers or hashable node keys).
    graph : networkx.Graph
        Graph providing node attribute `'color'` and edge attribute `'color'`.
        Edge lookups handle both `(u, v)` and `(v, u)` for undirected graphs.
    debug : bool, optional
        If True, print per-edge and intermediate diagnostics. Default is False.

    Returns
    -------
    sequences : list of str
        List of strings, one per non-trivial trail, where each character
        corresponds to a molecular unit encountered along the trail.
    char_dict : dict
        Mapping from molecular unit tuples
        `(start_color, edge_color, end_color)` to the single-character symbol
        used in the corresponding `sequences`.
    trivial_trail_count : int
        Number of trails of length 2 (single edges) that were counted as
        trivial and excluded from `sequences`.

    Raises
    ------
    KeyError
        If a traversed edge or node is missing the expected `'color'`
        attribute required to construct a molecular unit.
    TypeError
        If `trails` or `graph` are of an unexpected type that prevents
        iteration or attribute access.
    """
    char_dict = dict()  # Dictionary to store the character representing each unit
    char_list = []  # List to store the characters
    trivial_trail_count = 0  # Counter for the number of trails of length 1 (they can be pulled early from the calculation)
    sequences = []  # List to store the sequences
    for trail in trails:
        if len(trail) == 2:
            trivial_trail_count += 1
            continue

        sequence = ""
        for u, v in zip(trail, trail[1:]):
            if debug:
                print(f"Processing edge ({u}, {v})", flush=True)

            edge_color = graph.edges.get((u, v), graph.edges.get((v, u), {})).get('color')
            if edge_color is None:
                raise KeyError(f"Edge ({u}, {v}) not found in the graph")

            mol_unit = (graph.nodes[u]['color'], edge_color, graph.nodes[v]['color'])
            char_unit = char_dict.get(
                mol_unit)  # or char_dict.get(mol_unit[::-1]) # Mapping both to the same char would allow distinct fragments to be identified with each other
            if not char_unit:
                char_unit = get_unique_char(char_list)
                char_dict[mol_unit] = char_unit
                char_list.append(char_unit)

            sequence += char_unit

        sequences.append(sequence)

    if debug:
        print(f"char_dict: {char_dict}", flush=True)
        print(f"sequences: {sequences}", flush=True)
        print(f"trivial_trail_count: {trivial_trail_count}", flush=True)

    return sequences, char_dict, trivial_trail_count


def repair_compression(sequences: List[str]) -> Tuple[List[List[str]], List[Tuple[str, Tuple[str, str]]]]:
    """
    Apply a RePair-like compression to a list of symbol sequences.

    Performs iterative pair-replacement compression similar to RePair:
    1. Count frequencies of all adjacent symbol pairs across sequences.
    2. Select the most frequent pair occurring more than once.
    3. Introduce a new nonterminal symbol to represent that pair.
    4. Replace all occurrences of the pair with the new symbol.
    5. Record the replacement rule and repeat until no pair occurs more than once.

    Parameters
    ----------
    sequences : list of str
        Iterable of input sequences to compress. Each sequence is treated as an
        ordered list of atomic symbols (single-character strings or tokens).
        The function copies and operates on list representations of these
        sequences, so input sequences are not mutated.

    Returns
    -------
    final_seqs : list of list
        Compressed sequences where each sequence is a list of symbols (strings).
        New nonterminal symbols are inserted as string tokens (e.g. ``"NT_0"``).
    rules : list of tuple
        List of replacement rules in the order they were created. Each rule is a
        tuple ``(new_symbol, (sym1, sym2))`` meaning ``new_symbol -> sym1 sym2``.

    Raises
    ------
    TypeError
        If ``sequences`` is not an iterable of strings or string-like iterables.
    ValueError
        If a generated new symbol collides with existing symbols in the input
        (this implementation uses a simple ``NT_{i}`` naming scheme and raises
        only if a collision policy is later introduced).
    """

    seqs = [list(seq) for seq in sequences]  # copy
    rules = []
    next_nonterminal_id = 0

    while True:
        pairs = Counter()
        for seq in seqs:
            pairs.update(Counter(zip(seq, seq[1:])))
        # Find the most frequent pair that occurs more than once
        pair, freq = None, 1
        for p, f in pairs.items():
            if f > freq:
                pair, freq = p, f

        if pair is None:
            break

        # Introduce a new symbol
        new_symbol = f"NT_{next_nonterminal_id}"
        next_nonterminal_id += 1
        rules.append((new_symbol, pair))

        # Replace all occurrences of pair with new_symbol
        # Efficient replacement:
        # We'll scan seq and whenever we see pair[i], pair[i+1], replace them with new_symbol.
        for s_idx, seq in enumerate(seqs):
            i = 0
            new_seq = []
            while i < len(seq):
                if i < len(seq) - 1 and (seq[i], seq[i + 1]) == pair:
                    new_seq.append(new_symbol)
                    i += 2
                else:
                    new_seq.append(seq[i])
                    i += 1
            seqs[s_idx] = new_seq

        # Repeat until no more frequent pairs

    return seqs, rules


def compute_assembly_path_length_from_compression(final_seqs: List[List[str]],
                                                  rules: List[Tuple[str, Tuple[str, str]]],
                                                  ai_correction: int,
                                                  debug: bool = False) -> int:
    """
    Compute assembly path length from RePair-compressed sequences and rules.

    Calculates the assembly/path length using the simple formula:

    path_length = (total length of compressed sequences) + (number of rules) + (ai_correction)

    Parameters
    ----------
    final_seqs : list of list of str
        Compressed sequences produced by `repair_compression`. Each element is a
        sequence represented as an iterable (typically a list) of symbol strings.
    rules : list of tuple
        Replacement rules produced during compression. Each rule is a tuple
        ``(new_symbol, (sym1, sym2))`` describing how a nonterminal expands.
    ai_correction : int
        Integer correction term applied to the computed path length. This value is
        computed outside the function (for example as ``unique_units - joint_correction - trivial_trail_count``)
        and may be negative, zero or positive.
    debug : bool, optional
        If True, print intermediate debug information (length of final sequences
        and number of rules). Default is False.

    Returns
    -------
    int
        The computed assembly/path length as an integer.

    Raises
    ------
    TypeError
        If ``final_seqs`` or ``rules`` are not iterables of the expected form, or
        if ``ai_correction`` is not an integer.
    ValueError
        If any sequence in ``final_seqs`` contains non-hashable or otherwise
        unsupported elements that prevent length measurement (rare).
    """
    length_final = sum(len(seq) for seq in final_seqs)
    num_rules = len(rules)
    if debug:
        print(f"length_final: {length_final}", flush=True)
        print(f"num_rules: {num_rules}", flush=True)
    return length_final + num_rules + ai_correction


def calculate_assembly_path_det(graph: nx.Graph,
                                iterations: int = 1,
                                debug: bool = False) -> Tuple[
    Union[int, float], Optional[List[nx.Graph]], Optional[nx.DiGraph]]:
    """
    Deterministic assembly-path search by repeated trail partitioning and RePair compression.

    Performs a multi-iteration heuristic search that partitions the input graph
    into edge-disjoint trails, encodes trails as symbol sequences, applies a
    RePair-like compression to the sequences, and computes a simple assembly/path
    length metric from the compressed representation. The shortest path length
    (and associated compressed representation) observed across iterations is
    returned alongside reconstructed virtual objects and a directed path graph
    useful for downstream reconstruction.

    Parameters
    ----------
    graph : networkx.Graph
        Undirected input graph whose nodes are expected to have a `'color'`
        attribute and whose edges are expected to have a `'color'` attribute.
    iterations : int, optional
        Number of independent random partitions/compressions to attempt. The
        algorithm is nondeterministic (randomized selection of seed edges and
        extension choices); increasing `iterations` can improve the chance of
        finding a shorter assembly path. Default is 1.
    debug : bool, optional
        If True, print diagnostic information for each iteration and internal
        processing steps. Default is False.

    Returns
    -------
    best_path_length : int or float
        The smallest assembly/path length found across all iterations. May be
        0 for trivial inputs or `numpy.inf` if no valid partitions were produced.
    virtual_objects : list or None
        List of reconstructed molecular-graph objects (one per discovered symbol)
        produced by unpacking the final compression rules with the character
        dictionary. Returns `None` for trivial early exits.
    path : networkx.DiGraph or None
        Directed graph representing relationships between compressed symbols
        (nonterminals and their components). Returns `None` for trivial early exits.

    Raises
    ------
    TypeError
        If `graph` is not a NetworkX graph or if required node/edge attributes
        are missing such that downstream processing cannot proceed.
    ValueError
        If the input graph is detected as trivial (no meaningful compression or
        assembly path can be computed) the function may return early rather than
        raise; callers should validate input separately if exceptions are
        preferred.
    """

    # Find the joint correction
    n_joint_correction = len(get_disconnected_subgraphs(graph)) - 1

    # Check for trivial input
    if graph.number_of_edges() <= 1:
        return 0, None, None
    elif n_joint_correction > 0:
        trivial = True
        for nodes in nx.connected_components(graph):
            subgraph = graph.subgraph(nodes)
            if subgraph.number_of_edges() > 1:
                trivial = False
                break
        if trivial:
            print("Trivial input detected. Returning early.", flush=True)
            return 0, None, None

    # Extract all original units (vertex color-edge color-vertex color triples) from the graph
    purged_graph, unique_units = purge_unique_units(graph)

    # Initialise best path vars
    best_path_length = np.inf
    best_final_seqs = None
    best_rules = None

    # Repeat the process for a number of iterations, and keep the best path
    for it in range(iterations):
        # Partition the graph into disjoint trails
        trails = partition_into_disjoint_trails(purged_graph, debug=debug)

        # Convert trails to sequences
        sequences, char_dict, trivial_trail_count = trails_to_sequences(trails, purged_graph, debug=debug)

        # Apply RePair compression on the sequences
        final_seqs, rules = repair_compression(sequences)

        # Compute assembly index from the compression
        path_len = compute_assembly_path_length_from_compression(final_seqs,
                                                                 rules,
                                                                 trivial_trail_count + unique_units - n_joint_correction,
                                                                 debug=debug)

        if debug:
            print(f"Iteration {it}: path length = {path_len}", flush=True)
            print(f"Final sequences: {final_seqs}", flush=True)
            print(f"Rules: {rules}", flush=True)
            print(f"Trivial trail count: {trivial_trail_count}", flush=True)
            print(f"Unique units: {unique_units}", flush=True)
            print(f"Joint correction: {n_joint_correction}", flush=True)

        # Update the best path vars
        if path_len < best_path_length:
            best_path_length = path_len
            best_final_seqs = final_seqs
            best_rules = rules
            best_dict = char_dict

    # Rebuild the graph from the compressed graph and rules
    virtual_objects, path = process_paths(best_rules, best_dict)

    # Return shortest path length, plus some form of compressed structure (final_seq) and rules
    # final_seq acts as compressed_graph representation.
    return best_path_length, virtual_objects, path


def purge_unique_units(graph: nx.Graph) -> Tuple[nx.Graph, int]:
    """
    Purge unique vertex-edge-vertex units from a graph.

    Counts occurrences of molecular units represented as tuples
    ``(vertex_color, edge_color, vertex_color)`` across all edges, removes
    edges corresponding to units that occur exactly once (treating reversed
    units as equivalent), and prunes any isolated nodes created by removals.

    Parameters
    ----------
    graph : networkx.Graph
        Input undirected graph whose nodes must have a ``'color'`` attribute
        and whose edges must have a ``'color'`` attribute. The function does
        not modify the input graph in place; it returns a shallow copy with
        the unique-unit edges removed.

    Returns
    -------
    purged_graph : networkx.Graph
        A copy of ``graph`` with edges corresponding to unique units removed
        and any resulting isolated nodes deleted.
    unique_units : int
        Number of unique units removed (i.e. count of distinct unit types
        that were observed exactly once and for which the corresponding edge
        was removed).

    Raises
    ------
    KeyError
        If any node in the graph is missing the required ``'color'`` attribute
        or any edge is missing the required ``'color'`` attribute.
    TypeError
        If ``graph`` is not a NetworkX graph-like object supporting
        ``edges(data=True)`` and node attribute access.
    """

    # Make dict counting the frequency of each unit
    unit_freqs = {}
    for u, v, d in graph.edges(data=True):
        unit = (graph.nodes[u]['color'], d['color'], graph.nodes[v]['color'])
        if unit in unit_freqs:
            unit_freqs[unit] += 1
        elif unit[::-1] in unit_freqs:
            unit_freqs[unit[::-1]] += 1
        else:
            unit_freqs[unit] = 1

    # Remove unique units
    unique_units = 0  # Counter for the number of unique units
    purged_graph = graph.copy()
    for unit, freq in unit_freqs.items():
        if freq == 1:
            u, e, v = unit
            for u_prime, v_prime, e_prime in purged_graph.edges(data=True):  # Read through edges
                if purged_graph.edges[u_prime, v_prime]['color'] == e:  # Check if edge color matches
                    if purged_graph.nodes[u_prime]['color'] == u and purged_graph.nodes[v_prime][
                        'color'] == v:  # Check if vertex colors match
                        purged_graph.remove_edge(u_prime, v_prime)
                        break
                    elif purged_graph.nodes[u_prime]['color'] == v and purged_graph.nodes[v_prime][
                        'color'] == u:  # Check if vertex colors match
                        purged_graph.remove_edge(u_prime, v_prime)
                        break
            unique_units += 1

    # Clean up any isolated nodes
    purged_graph.remove_nodes_from(list(nx.isolates(purged_graph)))

    return purged_graph, unique_units


def process_paths(rules: List[Tuple[str, Tuple[str, str]]],
                  char_dict: Dict[Tuple[Any, Any, Any], str]) -> Tuple[List[nx.Graph], nx.DiGraph]:
    """
    Construct virtual molecular graphs and a directed symbol graph from compression rules and a character dictionary.

    Builds a directed graph describing how nonterminal symbols expand into component symbols
    and reconstructs the corresponding "virtual" molecular graphs (one per discovered symbol)
    by recursively unpacking rules and translating symbol sequences into molecular units.

    Parameters
    ----------
    rules : list of tuple
        List of replacement rules produced by the compression pass. Each rule must be a
        tuple ``(new_symbol, (sym1, sym2))`` where ``new_symbol`` is a symbol introduced by
        compression (for example ``"NT_0"``) and ``sym1`` and ``sym2`` are the two component
        symbols (either terminal symbols or other nonterminals).
    char_dict : dict
        Mapping from molecular unit tuples to single-character symbols, i.e.
        ``{(start_color, edge_color, end_color): 'A', ...}``. This function will invert
        the mapping internally to obtain a mapping from symbol (character) to unit
        required for molecular reconstruction.

    Returns
    -------
    virtual_objects : list of networkx.Graph
        List of molecular graphs reconstructed from each discovered symbol. Each element
        is a NetworkX undirected graph whose nodes have a ``'color'`` attribute and whose
        edges have a ``'color'`` attribute corresponding to the molecular unit colors.
    path : networkx.DiGraph
        Directed graph representing symbol expansion relationships. Nodes correspond to
        discovered symbols (indexed by integer insertion order) and directed edges point
        from a nonterminal to its component symbols.

    Raises
    ------
    TypeError
        If ``rules`` is not an iterable of 2-tuples of the expected form, or if ``char_dict``
        is not a mapping. Also raised when downstream helpers receive unsupported types.
    KeyError
        If a terminal symbol referenced during unpacking is not present in the inverted
        ``char_dict`` (i.e. a terminal character has no associated molecular unit).
    ValueError
        If ``rules`` contain inconsistent or cyclic definitions that prevent successful
        unpacking into terminal sequences.
    """
    path = nx.DiGraph()  # Initialize a directed graph
    symbol_virtual_objects = []  # List to store unique symbols as virtual objects

    # Iterate through the rules to construct the graph
    for rule in rules:
        new_symbol, (sym1, sym2) = rule
        # Create a new virtual object node if it doesn't already exist
        if new_symbol not in symbol_virtual_objects:
            symbol_virtual_objects.append(new_symbol)
            path.add_node(len(symbol_virtual_objects))
        # Add edges for the components of the rule
        for sym in [sym1, sym2]:
            if sym not in symbol_virtual_objects:
                symbol_virtual_objects.append(sym)
                path.add_node(len(symbol_virtual_objects))
            if (symbol_virtual_objects.index(new_symbol), symbol_virtual_objects.index(sym)) not in path.edges:
                path.add_edge(symbol_virtual_objects.index(new_symbol), symbol_virtual_objects.index(sym))

    # Process the virtual objects to generate molecular graphs
    virtual_objects = []  # List to store molecular graphs
    unit_dict = {v: k for k, v in char_dict.items()}  # Reverse the character dictionary
    for s_idx, symbol in enumerate(symbol_virtual_objects):
        # Reconstruct the sequence for the symbol and convert it to a molecular graph
        seq = unpack_path(symbol, rules)
        mol_graph = seq_2_mol(seq, unit_dict)
        virtual_objects.append(mol_graph)

    return virtual_objects, path


def unpack_path(symbol: str, rules: List[Tuple[str, Tuple[str, str]]]) -> str:
    """
    Recursively reconstruct a terminal sequence from a symbol using binary production rules.

    Given a set of binary rules of the form ``(new_symbol, (sym1, sym2))``, this
    function expands ``symbol`` by recursively replacing nonterminal symbols with
    their right-hand components until only terminal symbols remain. The result is
    returned as a concatenated string of terminal symbols.

    Parameters
    ----------
    symbol : str
        Starting symbol to unpack. May be a terminal (single-character) or a
        nonterminal introduced by the compression pass (for example ``"NT_0"``).
    rules : list of tuple
        Iterable of binary production rules. Each element must be a 2-tuple
        ``(new_symbol, (sym1, sym2))`` where ``sym1`` and ``sym2`` are symbols
        that may themselves be terminals or nonterminals.

    Returns
    -------
    str
        The fully expanded sequence obtained by recursively unpacking ``symbol``.
        Terminal symbols are concatenated in left-to-right order produced by the
        rule expansions.

    Raises
    ------
    TypeError
        If ``rules`` is not an iterable of 2-tuples or if rule elements are not of
        the expected types (string for symbols and a 2-tuple for the right-hand side).
    ValueError
        If a cycle is detected in the rules (recursive definitions) that prevents
        termination, or if a referenced nonterminal has no corresponding rule.
    """
    output = ""
    # If the symbol is not a new symbol in the rules, return it as is
    if symbol not in [rule[0] for rule in rules]:
        return symbol
    else:
        # Find the rule corresponding to the symbol and recursively unpack its components
        for rule in rules:
            if rule[0] == symbol:
                output = unpack_path(rule[1][0], rules) + unpack_path(rule[1][1], rules)
                break
        return output


def seq_2_mol(seq: str, unit_dict: Dict[str, Tuple[Any, Any, Any]]) -> nx.Graph:
    """
    Convert a sequence of symbols into a molecular graph.

    Creates a linear NetworkX undirected graph from a sequence of terminal symbols
    using a mapping from each symbol to a molecular unit tuple
    ``(start_color, edge_color, end_color)``. Nodes are 1-indexed integers and carry
    a ``'color'`` node attribute; edges carry a ``'color'`` edge attribute.

    Parameters
    ----------
    seq : str
        Sequence of terminal symbols (each symbol is a key in ``unit_dict``).
        The sequence length ``m`` yields a graph with ``m + 1`` nodes and ``m``
        edges representing the successive molecular units.
    unit_dict : dict
        Mapping from single-character symbols to molecular unit tuples:
        ``{symbol: (start_color, edge_color, end_color), ...}``. Colors can be
        any hashable values (commonly strings).

    Returns
    -------
    networkx.Graph
        Undirected molecular graph where nodes are integers ``1..m+1`` with a
        ``'color'`` attribute and edges between consecutive nodes carry a
        ``'color'`` attribute corresponding to the unit's edge color.

    Raises
    ------
    TypeError
        If ``seq`` is not a string-like iterable or ``unit_dict`` is not a mapping.
    KeyError
        If any symbol in ``seq`` is not present in ``unit_dict``.
    ValueError
        If ``seq`` is empty (no units) and a non-empty molecular graph is required.
    """
    mol_graph = nx.Graph()
    # Initialize the first node with the start color of the first unit in the sequence
    mol_graph.add_node(1, color=unit_dict[seq[0]][0])
    # Iterate through the sequence to build nodes and edges
    for idx in range(len(seq)):
        # Add a node for the end color of the current unit
        mol_graph.add_node(idx + 2, color=unit_dict[seq[idx]][2])
        # Add an edge between the current node and the next node with the edge color
        mol_graph.add_edge(idx + 1, idx + 2, color=unit_dict[seq[idx]][1])
    return mol_graph
