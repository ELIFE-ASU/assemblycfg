import collections
import string

import networkx as nx


def rules_to_graph(rules, virt_obj):
    """
    Converts a list of rules into a directed graph and adds virtual objects as nodes.

    Args:
        rules (list): Rules in the format "A + B = C".
        virt_obj (list): Virtual objects to be added as nodes.

    Returns:
        networkx.DiGraph: Directed graph with rules and virtual objects.
    """
    # Create a directed graph and add virtual objects as nodes
    graph = nx.DiGraph()
    graph.add_nodes_from(virt_obj)

    # Add edges based on rules
    for rule in rules:
        a, b, c = rule.replace(" + ", " = ").split(" = ")
        graph.add_edge(a, c)
        graph.add_edge(b, c)

    return graph


def repair(s):
    """
    Repairs the input string by replacing the most frequent adjacent pairs of symbols with new non-terminal symbols.

    Args:
        s (str): The input string to be repaired.

    Returns:
        tuple: A tuple containing:
            - symbols (list): The list of symbols after replacement.
            - productions (dict): A dictionary of productions where keys are new non-terminal symbols and values are the pairs they replace.
    """
    symbols = list(s)
    productions = {}
    non_terminal_counter = 1

    while True:
        # Count the frequency of adjacent pairs and filter those occurring more than once
        pair_counts = collections.Counter(zip(symbols, symbols[1:]))
        frequent_pairs = {pair: count for pair, count in pair_counts.items() if count > 1}

        if not frequent_pairs:
            break

        # Find the most frequent pair and create a new non-terminal
        most_frequent_pair = max(frequent_pairs, key=frequent_pairs.get)
        new_non_terminal = f'A{non_terminal_counter}'
        non_terminal_counter += 1
        productions[new_non_terminal] = list(most_frequent_pair)

        i = 0
        while i < len(symbols) - 1:
            # Check if the current pair matches the most frequent pair
            if (symbols[i], symbols[i + 1]) == most_frequent_pair:
                # Replace the pair with the new non-terminal
                symbols[i:i + 2] = [new_non_terminal]
                i = max(i - 1, 0)  # Step back to handle overlapping pairs
            else:
                i += 1

    return symbols, productions


def convert_to_cnf(start_symbol, productions):
    """
    Converts a context-free grammar (CFG) to Chomsky Normal Form (CNF).

    Args:
        start_symbol (str): The start symbol of the CFG.
        productions (dict): A dictionary where keys are non-terminal symbols and values are lists of symbols.

    Returns:
        tuple: The new start symbol and a dictionary of CNF productions.
    """
    cnf_productions = {}
    new_nt_counter = 1

    # Map terminals to new non-terminals
    terminals = {s for exp in productions.values() for s in exp if s in string.ascii_lowercase}
    terminals.update(s for s in start_symbol if s in string.ascii_lowercase)
    terminal_map = {t: f'T_{t}' for t in terminals}
    cnf_productions.update({nt: [t] for t, nt in terminal_map.items()})

    def replace_terminals(symbols):
        return [terminal_map.get(s, s) for s in symbols]

    # Convert productions to CNF
    for nt, expansion in productions.items():
        expansion = replace_terminals(expansion)
        while len(expansion) > 2:
            new_nt = f'N{new_nt_counter}'
            cnf_productions[new_nt] = expansion[:2]
            expansion = [new_nt] + expansion[2:]
            new_nt_counter += 1
        cnf_productions[nt] = expansion

    # Handle the start symbol
    start_symbols = replace_terminals(start_symbol)
    start_nt = 'S'
    while len(start_symbols) > 2:
        new_nt = f'N{new_nt_counter}'
        cnf_productions[new_nt] = start_symbols[:2]
        start_symbols = [new_nt] + start_symbols[2:]
        new_nt_counter += 1
    cnf_productions[start_nt] = start_symbols

    return start_nt, cnf_productions


def ai_core(s):
    """
    Converts the input string into Chomsky Normal Form (CNF) and calculates the production count.

    Args:
        s (str): The input string to be processed.

    Returns:
        tuple: A tuple containing:
            - production_count (int): The number of productions in the CNF minus the number of unique symbols in the input string.
            - cnf_productions (dict): A dictionary of CNF productions.
    """
    start_symbol, productions = repair(s)
    start_nt, cnf_productions = convert_to_cnf(start_symbol, productions)
    return len(cnf_productions) - len(set(s)), cnf_productions


def get_rules(s, production, f_print=False):
    """
    Generates a list of rules from the given production dictionary by performing a topological sort.

    Args:
        s (str): The input string to be processed.
        production (dict): A dictionary where keys are non-terminal symbols and values are lists of symbols.
        f_print (bool): Flag to print the rules and path length.

    Returns:
        list: A list of rules in the format "A + B = C".
    """
    in_degrees = collections.defaultdict(int)
    adj = collections.defaultdict(list)
    tmap = {}

    # Initialize in-degrees and adjacency list
    for c in s:
        in_degrees[c] = 0
    for course, prereqs in production.items():
        for req in prereqs:
            in_degrees[course] += 1
            adj[req].append(course)

    # Start queue with symbols having zero in-degrees
    start_q = collections.deque(symbol for symbol, ins in in_degrees.items() if ins == 0)
    if f_print:
        print(f"Processing {s}", flush=True)
        print(f"Start symbols: {', '.join(start_q)}", flush=True)
        print("Joins:", flush=True)

    # Perform topological sort
    rules = []
    while start_q:
        symbol = start_q.popleft()
        tmap[symbol] = symbol
        for neighbor in adj[symbol]:
            in_degrees[neighbor] -= 1
            if in_degrees[neighbor] == 0:
                start_q.append(neighbor)

        if symbol in production:
            if len(production[symbol]) == 2:
                a, b = production[symbol]
                tmap[symbol] = tmap[a] + tmap[b]
                rules.append(f"{tmap[a]} + {tmap[b]} = {tmap[symbol]}")
            else:
                tmap[symbol] = production[symbol][0]

    if f_print:
        print("\n".join(rules), flush=True)

    return rules


def extract_virtual_objects(rules):
    """
    Extracts a sorted set of objects from the given list of rules, excluding the entry on the right side of the last "=".

    Args:
        rules (list): A list of rules in the format "A + B = C".

    Returns:
        list: A list of objects found in the equations, sorted by the size of the string, excluding the last result.
    """
    if not rules:
        return []

    # Extract all objects from the rules
    objects = {obj for rule in rules for obj in rule.replace(" + ", " = ").split(" = ")}

    # Remove the last result
    last_result = rules[-1].split(" = ")[-1]
    objects.discard(last_result)

    # Return the sorted list of objects by string length
    return sorted(objects, key=len)


def ai_with_pathways(s, f_print=False):
    """
    Takes the production rules from ai_upper. Performs a topological sort to find order
    of join operations.

    Args:
        s (str): input string to be processed.
        f_print (bool): flag to print the rules and path length.

    Returns:
        ai_count (int): the final path length.
        virt_obj (list): a list of virtual objects found in the equations, sorted by the size of the string.
        rules (list): a list of rules in the format "A + B = C".
    """
    # Get the production rules and path length
    ai_count, production = ai_core(s)
    # Get the rules
    rules = get_rules(s, production, f_print=f_print)
    # Extract virtual objects
    virt_obj = extract_virtual_objects(rules)
    return ai_count, virt_obj, rules_to_graph(rules, virt_obj)
