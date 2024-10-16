import collections

def repair(s):
    symbols = list(s)
    productions = {}
    non_terminal_counter = 1

    while True:
        pair_counts = collections.Counter()
        for i in range(len(symbols) - 1):
            pair = (symbols[i], symbols[i + 1])
            pair_counts[pair] += 1

        frequent_pairs = {pair: count for pair, count in pair_counts.items() if count > 1}

        if not frequent_pairs:
            break

        most_frequent_pair = max(frequent_pairs, key=frequent_pairs.get)

        new_non_terminal = f'A{non_terminal_counter}'
        non_terminal_counter += 1

        productions[new_non_terminal] = list(most_frequent_pair)

        i = 0
        while i < len(symbols) - 1:
            if symbols[i] == most_frequent_pair[0] and symbols[i + 1] == most_frequent_pair[1]:
                symbols[i] = new_non_terminal
                del symbols[i + 1]
                if i > 0:
                    i -= 1
            else:
                i += 1

    start_symbol = symbols
    return start_symbol, productions

def convert_to_cnf(start_symbol, productions):
    cnf_productions = {}
    new_non_terminal_counter = 1

    # Step 1: Create non-terminals for terminals
    terminals = set()
    for expansion in productions.values():
        for symbol in expansion:
            if len(symbol) == 1 and symbol.islower():
                terminals.add(symbol)
    for symbol in start_symbol:
        if len(symbol) == 1 and symbol.islower():
            terminals.add(symbol)

    terminal_non_terminals = {t: f'T_{t}' for t in terminals}

    # Add terminal productions
    for t, nt in terminal_non_terminals.items():
        cnf_productions[nt] = [t]

    # Step 2: Replace terminals in productions
    def replace_terminals(symbols):
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

def AI_UPPER(s):
    start_symbol, productions = repair(s)
    start_nt, cnf_productions = convert_to_cnf(start_symbol, productions)
    production_count = len(cnf_productions) - len(set(s))
    return production_count


s = "abracadabra"
print(AI_UPPER(s))
