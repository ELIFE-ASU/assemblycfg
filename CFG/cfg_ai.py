import collections


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


def ai_upper(s):
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


def ai_upper_with_pathways(s):
    """
    Takes the production rules from ai_upper. Performs a topological sort to find order
    of join operations. Prints the proper AI joins in order.

    Args:
        s (str): input string to be processed.

    Returns:
        ai_count (int): the final path length.
    """
    print(f"Processing {s}", flush=True)
    ai_count, production = ai_upper(s)
    in_degrees = collections.defaultdict(int)
    adj = collections.defaultdict(list)
    for c in s:
        in_degrees[c] = 0
        adj[c] = []

    tmap = collections.defaultdict(str)

    for course, prereqs in production.items():
        for req in prereqs:
            in_degrees[course] += 1
            adj[req].append(course)

    start_q = collections.deque()
    for symbol, ins in in_degrees.items():
        if ins == 0:
            start_q.append(symbol)

    print(f"START SYMBOLS: {','.join(start_q)}", flush=True)
    print("JOINS: ", flush=True)
    rules = []
    q = collections.deque()
    while start_q:
        symbol = start_q.popleft()
        tmap[symbol] = symbol
        for neighbor in adj[symbol]:
            in_degrees[neighbor] -= 1
            if in_degrees[neighbor] == 0:
                q.append(neighbor)
    while q:
        symbol = q.popleft()
        if len(production[symbol]) == 2:
            a, b = production[symbol]
            tmap[symbol] = tmap[a] + tmap[b]
            rule = f"{tmap[a]} + {tmap[b]} = {tmap[symbol]}"
            rules.append(rule)
            print(rule, flush=True)
        else:
            tmap[symbol] = production[symbol][0]
        for neighbor in adj[symbol]:
            in_degrees[neighbor] -= 1
            if in_degrees[neighbor] == 0:
                q.append(neighbor)
    print(f"PATH LENGTH: {ai_count}", flush=True)
    return ai_count, rules
