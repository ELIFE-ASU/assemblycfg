import collections

import networkx as nx


def rules_to_graph(rules, virt_obj):
    """
    Converts a list of rules into a directed NetworkX graph and adds virtual objects as nodes.

    Args:
        rules (list): A list of rules in the format "A + B = C".
        virt_obj (list): A list of virtual objects to be added as nodes in the graph.

    Returns:
        networkx.DiGraph: A directed graph representing the rules with virtual objects as nodes.
    """
    # Split each rule into parts by replacing " + " with " = " and then splitting by " = "
    parts = [rule.replace(" + ", " = ").split(" = ") for rule in rules]

    # Create a directed graph
    G = nx.DiGraph()
    # Add nodes for each virtual object
    G.add_nodes_from(virt_obj)

    # Add edges based on parts
    for part in parts:
        G.add_edge(part[0], part[2])
        G.add_edge(part[1], part[2])

    return G


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
        # Count the frequency of each adjacent pair of symbols
        pair_counts = collections.Counter((symbols[i], symbols[i + 1]) for i in range(len(symbols) - 1))
        # Filter pairs that occur more than once
        frequent_pairs = {pair: count for pair, count in pair_counts.items() if count > 1}

        if not frequent_pairs:
            break

        # Find the most frequent pair
        most_frequent_pair = max(frequent_pairs, key=frequent_pairs.get)
        new_non_terminal = f'A{non_terminal_counter}'
        non_terminal_counter += 1
        # Add the new non-terminal to the productions dictionary
        productions[new_non_terminal] = list(most_frequent_pair)

        i = 0
        while i < len(symbols) - 1:
            # Replace occurrences of the most frequent pair with the new non-terminal
            if symbols[i:i + 2] == list(most_frequent_pair):
                symbols[i] = new_non_terminal
                del symbols[i + 1]
                i = max(i - 1, 0)
            else:
                i += 1

    return symbols, productions


def convert_to_cnf(start_symbol, productions):
    """
    Converts a context-free grammar (CFG) to Chomsky Normal Form (CNF).

    Args:
        start_symbol (str): The start symbol of the CFG.
        productions (dict): A dictionary where keys are non-terminal symbols and values are lists of symbols (terminals or non-terminals).

    Returns:
        tuple: A tuple containing:
            - start_nt (str): The new start symbol for the CNF.
            - cnf_productions (dict): A dictionary of CNF productions.
    """
    cnf_productions = {}
    new_non_terminal_counter = 1

    # Step 1: Create non-terminals for terminals
    terminals = {symbol for expansion in productions.values() for symbol in expansion if
                 len(symbol) == 1 and symbol.islower()}
    terminals.update(symbol for symbol in start_symbol if len(symbol) == 1 and symbol.islower())

    terminal_non_terminals = {t: f'T_{t}' for t in terminals}

    # Add terminal productions
    cnf_productions.update({nt: [t] for t, nt in terminal_non_terminals.items()})

    # Step 2: Replace terminals in productions
    def replace_terminals(symbols):
        """
        Replaces terminal symbols in a list of symbols with their corresponding non-terminals.

        Args:
            symbols (list): A list of symbols (terminals or non-terminals).

        Returns:
            list: A list of symbols with terminals replaced by their corresponding non-terminals.
        """
        return [terminal_non_terminals.get(s, s) for s in symbols]

    # Process existing productions
    for nt, expansion in productions.items():
        new_expansion = replace_terminals(expansion)
        while len(new_expansion) > 2:
            new_nt = f'N{new_non_terminal_counter}'
            new_non_terminal_counter += 1
            cnf_productions[new_nt] = new_expansion[:2]
            new_expansion = [new_nt] + new_expansion[2:]
        cnf_productions[nt] = new_expansion

    # Handle the start symbol
    start_symbols = replace_terminals(start_symbol)
    start_nt = 'S'
    while len(start_symbols) > 2:
        new_nt = f'N{new_non_terminal_counter}'
        new_non_terminal_counter += 1
        cnf_productions[new_nt] = start_symbols[:2]
        start_symbols = [new_nt] + start_symbols[2:]
    cnf_productions[start_nt] = start_symbols

    return start_nt, cnf_productions


def ai_core(s):
    """
    Processes the input string to convert it into Chomsky Normal Form (CNF) and calculates the production count.

    Args:
        s (str): The input string to be processed.

    Returns:
        tuple: A tuple containing:
            - production_count (int): The number of productions in the CNF minus the number of unique symbols in the input string.
            - cnf_productions (dict): A dictionary of CNF productions.
    """
    start_symbol, productions = repair(s)
    start_nt, cnf_productions = convert_to_cnf(start_symbol, productions)
    production_count = len(cnf_productions) - len(set(s))
    return production_count, cnf_productions


def get_rules(s, production, f_print=False):
    """
    Generates a list of rules from the given production dictionary by performing a topological sort.

    Args:
        s (str): The input string to be processed.
        production (dict): A dictionary where keys are non-terminal symbols and values are lists of symbols (terminals or non-terminals).
        f_print (bool): Flag to print the rules and path length.

    Returns:
        list: A list of rules in the format "A + B = C".
    """
    in_degrees = collections.defaultdict(int)
    adj = collections.defaultdict(list)
    tmap = collections.defaultdict(str)

    # Initialize in-degrees for each symbol in the input string
    for c in s:
        in_degrees[c] = 0

    # Build the adjacency list and in-degrees for each course
    for course, prereqs in production.items():
        for req in prereqs:
            in_degrees[course] += 1
            adj[req].append(course)

    # Initialize the start queue with symbols having zero in-degrees
    start_q = collections.deque([symbol for symbol, ins in in_degrees.items() if ins == 0])
    if f_print:
        print(f"Processing {s}", flush=True)
        print(f"START SYMBOLS: {','.join(start_q)}", flush=True)
        print("JOINS: ", flush=True)

    q = collections.deque()
    while start_q:
        symbol = start_q.popleft()
        tmap[symbol] = symbol
        for neighbor in adj[symbol]:
            in_degrees[neighbor] -= 1
            if in_degrees[neighbor] == 0:
                q.append(neighbor)

    rules = []
    while q:
        symbol = q.popleft()
        if len(production[symbol]) == 2:
            a, b = production[symbol]
            tmap[symbol] = tmap[a] + tmap[b]
            rules.append(f"{tmap[a]} + {tmap[b]} = {tmap[symbol]}")
        else:
            tmap[symbol] = production[symbol][0]
        for neighbor in adj[symbol]:
            in_degrees[neighbor] -= 1
            if in_degrees[neighbor] == 0:
                q.append(neighbor)

    if f_print:
        for rule in rules:
            print(rule, flush=True)

    return rules


def extract_virtual_objects(rules):
    """
    Extracts a sorted set of objects from the given list of rules, excluding the entry on the right side of the last "=".

    Args:
        rules (list): A list of rules in the format "A + B = C".

    Returns:
        list: A list of objects found in the equations, sorted by the size of the string, excluding the last result.
    """
    objects = set()
    for rule in rules:
        # Replace " + " with " = " and split the rule into individual objects
        parts = rule.replace(" + ", " = ").split(" = ")
        # Add the objects to the set
        objects.update(parts)
    # Remove the last result
    if rules:
        last_result = rules[-1].split(" = ")[-1]
        objects.discard(last_result)
    # Return the sorted list of objects by the size of the string
    return sorted(objects, key=len)


def ai_with_pathways(s, f_print=False):
    """
    Takes the production rules from ai_upper. Performs a topological sort to find order
    of join operations. Prints the proper AI joins in order.

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
    if f_print:
        print(f"PATH LENGTH: {ai_count}", flush=True)
        print(f"VIRTUAL OBJECTS: {virt_obj}", flush=True)
    return ai_count, virt_obj, rules_to_graph(rules, virt_obj)
