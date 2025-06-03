import random
from collections import Counter
import networkx as nx
import numpy as np
from .utils import get_disconnected_subgraphs


def partition_into_disjoint_trails(graph, debug=False):
    """
    Partitions a NetworkX graph into disjoint trails (walks that do not reuse edges).

    Parameters:
        graph (networkx.Graph): The input graph.
        debug (bool): If True, prints debug information. Default is False.

    Returns:
        list: A list of disjoint trails, where each trail is a list of nodes.
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


def get_unique_char(input_str):
    """
    Find a unique character that is not present in the input string.

    Args:
        input_str (str or list): The input string or list of characters to check against.

    Returns:
        str: A unique character not present in the input string.
    """
    for i in range(1, 1114111):
        char = chr(i)
        if char not in input_str:
            return char


def trails_to_sequences(trails, graph, debug=False):
    """
    Convert a list of trails to a string with a corresponding assembly index to the parent graph.

    Parameters:
        trails (list): A list of trails, where each trail is a list of nodes.
        graph (networkx.Graph): The input graph.
        debug (bool): If True, prints debug information. Default is False.

    Returns:
        tuple: A tuple containing:
            - list: Strings representing the trails.
            - dict: A dictionary mapping the string units to the molecular units.
            - int: The number of extra edges (trails of length 1).
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


def repair_compression(sequences):
    """
    Apply a RePair-like compression on sequences (a list of strings).

    Steps:
    1. Count frequencies of all adjacent pairs.
    2. Choose the most frequent pair.
    3. Replace all occurrences of that pair with a new symbol.
    4. Record the rule for that pair.
    5. Repeat until no pair occurs more than once.

    Parameters:
        sequences (list of str): The list of sequences to be compressed.

    Returns:
        tuple: A tuple containing:
            - list of list of str: The compressed sequences.
            - list of tuple: The list of rules, each rule is (new_symbol, (sym1, sym2)).

    Symbols can be strings or tuples. We'll keep them as is. New symbols will be generated as special strings.
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


def compute_assembly_path_length_from_compression(final_seqs, rules, ai_correction, debug=False):
    """
    Compute assembly index using the formula:
    path_length = (length_of_final_sequence) + (number_of_rules) + (ai_correction)

    Parameters:
        final_seqs (list of list of str): The compressed sequences after RePair.
        rules (list of tuple): The list of rules generated during compression.
        ai_correction (int): The assembly index correction factor.
        debug (bool): If True, prints debug information. Default is False.

    Returns:
        int: The computed assembly path length.

    ai_correction is calculated as:
    number of unique units - joint ai correction - trivial trail count
    """
    length_final = sum(len(seq) for seq in final_seqs)
    num_rules = len(rules)
    if debug:
        print(f"length_final: {length_final}", flush=True)
        print(f"num_rules: {num_rules}", flush=True)
    return length_final + num_rules + ai_correction


def calculate_assembly_path_det(graph, iterations=1, debug=False):
    """
    Main function:
    1. Extract all original units (vertex color-edge color-vertex color triple) from the graph.
    2. Build a sequence from the graph (by partitioning graph into disjoint trails).
    3. Apply RePair compression on the sequence.
    4. Compute assembly index from the formula.
    5. Repeat from step 2 for a number of iterations and return the shortest path found.

    Also returns the compressed sequence and rules for reconstruction.

    Parameters:
        graph (networkx.Graph): The input graph.
        iterations (int): The number of iterations to perform. Default is 1.
        debug (bool): If True, prints debug information. Default is False.

    Returns:
        tuple: A tuple containing:
            - int: The shortest path length found.
            - dict: Virtual objects (currently empty).
            - list: The rules for reconstruction.
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


def purge_unique_units(graph):
    """
    Purge the unique units (vertex color-edge color-vertex color triples) from the graph.

    Parameters:
        graph (networkx.Graph): The input graph.

    Returns:
        tuple: A tuple containing:
            - networkx.Graph: The purged graph with unique units removed.
            - int: The number of unique units removed.
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


def process_paths(rules, char_dict):
    """
    Processes the rules and character dictionary to construct a directed graph and a list of molecular graphs.

    Parameters:
        rules (list of tuple): A list of rules, where each rule is a tuple in the form (new_symbol, (sym1, sym2)).
                               Each rule defines how a new symbol is composed of two other symbols.
        char_dict (dict): A dictionary mapping symbols to their corresponding molecular units.

    Returns:
        tuple: A tuple containing:
            - list: A list of molecular graphs (virtual objects) corresponding to the symbols.
            - networkx.DiGraph: A directed graph representing the relationships between symbols.
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


def unpack_path(symbol, rules):
    """
    Recursively reconstructs a sequence from a given symbol and a set of rules.

    Parameters:
        symbol (str): The starting symbol to unpack.
        rules (list of tuple): A list of rules, where each rule is a tuple in the form (new_symbol, (sym1, sym2)).
                               Each rule defines how a new symbol is composed of two other symbols.

    Returns:
        str: The reconstructed sequence as a string.
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


def seq_2_mol(seq, unit_dict):
    """
    Convert a sequence of symbols into a molecular graph.

    Parameters:
        seq (str): A string representing the sequence of symbols.
        unit_dict (dict): A dictionary mapping each symbol to its corresponding molecular unit.
                          Each molecular unit is a tuple (start_color, edge_color, end_color).

    Returns:
        networkx.Graph: A molecular graph where nodes represent atoms (with colors) and edges represent bonds (with colors).
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
