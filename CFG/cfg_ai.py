import collections


def repair(s):
    symbols = list(s)
    productions = {}
    non_terminal_counter = 1

    while True:
        pair_counts = collections.Counter((symbols[i], symbols[i + 1]) for i in range(len(symbols) - 1))
        frequent_pairs = {pair: count for pair, count in pair_counts.items() if count > 1}

        if not frequent_pairs:
            break

        most_frequent_pair = max(frequent_pairs, key=frequent_pairs.get)
        new_non_terminal = f'A{non_terminal_counter}'
        non_terminal_counter += 1
        productions[new_non_terminal] = list(most_frequent_pair)

        i = 0
        while i < len(symbols) - 1:
            if symbols[i:i + 2] == list(most_frequent_pair):
                symbols[i] = new_non_terminal
                del symbols[i + 1]
                i = max(i - 1, 0)
            else:
                i += 1

    return symbols, productions


def convert_to_cnf(start_symbol, productions):
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
    start_symbol, productions = repair(s)
    start_nt, cnf_productions = convert_to_cnf(start_symbol, productions)
    production_count = len(cnf_productions) - len(set(s))
    return production_count, cnf_productions
