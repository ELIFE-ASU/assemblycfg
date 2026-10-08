"""
Upper bounds on string assembly index from the RePair grammar.

RePair compresses the input into a straight-line grammar; converting it to
Chomsky normal form gives one join per binary production, which is an
assembly pathway.
"""
import collections
import heapq
import string
from typing import Callable, Hashable, List, Sequence, Tuple, Dict, Union

import networkx as nx


def repair_pairs(sequences: Sequence[Sequence[Hashable]],
                 new_symbol: Callable[[int], Hashable]
                 ) -> Tuple[List[List[Hashable]], List[Tuple[Hashable, Tuple[Hashable, Hashable]]]]:
    """
    Run RePair over symbol sequences in near-linear time.

    Repeatedly replaces the adjacent pair that occurs most often (counting
    overlapping occurrences, and never across sequences) with a fresh
    symbol, until no pair occurs more than once. Ties go to the pair that
    occurs first, and each replacement scans left to right. This is the
    textbook quadratic loop, made incremental: the sequences are linked
    lists, every pair keeps a heap of its positions, and the pairs sit in a
    heap keyed by count and first position, so each round touches only the
    occurrences it replaces and their neighbours.

    Parameters
    ----------
    sequences : sequence of sequences of hashable
        The input sequences. They are not modified.
    new_symbol : callable
        Maps the 0-based index of a new rule to its symbol, which must not
        occur in the input.

    Returns
    -------
    final_seqs : list of list
        The compressed sequences.
    rules : list of tuple
        ``(symbol, (left, right))`` for each rule, in creation order.
    """
    # Symbols are interned as ints: inputs first, then rules in order.
    names: List[Hashable] = []
    ids: Dict[Hashable, int] = {}
    sym: List[int] = []
    nxt: List[int] = []
    prv: List[int] = []
    starts: List[int] = []
    lengths: List[int] = []
    for seq in sequences:
        start = len(sym)
        starts.append(start)
        for symbol in seq:
            if symbol not in ids:
                ids[symbol] = len(names)
                names.append(symbol)
            sym.append(ids[symbol])
        end = len(sym)
        lengths.append(end - start)
        nxt.extend(range(start + 1, end + 1))
        prv.extend(range(start - 1, end - 1))
        if end > start:
            nxt[-1] = -1
            prv[start] = -1

    # Positions are never reordered, so a pair's first occurrence is its
    # smallest position. Each position holds at most one pair in its
    # lifetime, so position heaps need no deduplication.
    count: Dict[Tuple[int, int], int] = collections.defaultdict(int)
    where: Dict[Tuple[int, int], List[int]] = collections.defaultdict(list)
    for p, q in enumerate(nxt):
        if q != -1:
            pair = (sym[p], sym[q])
            count[pair] += 1
            where[pair].append(p)  # ascending, so already a heap

    def first(pair):
        heap = where[pair]
        while heap:
            p = heap[0]
            q = nxt[p]
            if q != -1 and sym[p] == pair[0] and sym[q] == pair[1]:
                return p
            heapq.heappop(heap)
        return -1

    queue = [(-c, where[pair][0], pair) for pair, c in count.items() if c > 1]
    heapq.heapify(queue)

    rules: List[Tuple[int, Tuple[int, int]]] = []
    while queue:
        c, p, pair = heapq.heappop(queue)
        if -c != count[pair] or p != first(pair):
            continue  # stale entry
        a, b = pair
        new = len(names)
        names.append(new_symbol(len(rules)))
        rules.append((new, pair))

        dirty = set()
        occurrences = where.pop(pair)
        occurrences.sort()
        for p in occurrences:
            q = nxt[p]
            # Skip occurrences consumed or altered earlier in this round.
            if q == -1 or sym[p] != a or sym[q] != b:
                continue
            left, right = prv[p], nxt[q]
            if left != -1:
                old = (sym[left], a)
                count[old] -= 1
                dirty.add(old)
            count[pair] -= 1
            if right != -1:
                old = (b, sym[right])
                count[old] -= 1
                dirty.add(old)
            # Merge q into p.
            sym[p] = new
            nxt[p] = right
            if right != -1:
                prv[right] = p
            nxt[q] = -1
            if left != -1:
                formed = (sym[left], new)
                count[formed] += 1
                heapq.heappush(where[formed], left)
                dirty.add(formed)
            if right != -1:
                formed = (new, sym[right])
                count[formed] += 1
                heapq.heappush(where[formed], p)
                dirty.add(formed)

        dirty.discard(pair)
        for other in dirty:
            if count[other] > 1:
                heapq.heappush(queue, (-count[other], first(other), other))

    final_seqs = []
    for start, length in zip(starts, lengths):
        # A sequence's first position is never merged away.
        seq = []
        p = start if length else -1
        while p != -1:
            seq.append(names[sym[p]])
            p = nxt[p]
        final_seqs.append(seq)
    return final_seqs, [(names[n], (names[x], names[y])) for n, (x, y) in rules]


def rules_to_graph(rules: List[str],
                   virt_obj: List[str]) -> nx.DiGraph:
    """
    Convert a list of join rules into a NetworkX directed graph and add virtual nodes.

    Parses rules in the form ``"A + B = C"`` and adds directed edges ``A -> C`` and
    ``B -> C`` for each rule. Virtual objects are added to the graph as isolated
    nodes (if not already present).

    Parameters
    ----------
    rules : list of str
        Sequence of rules where each rule is expected to contain two operands and
        a result using the literal separators ``' + '`` and `` = '`` (e.g.
        ``"A + B = C"``). Empty list is accepted and yields a graph containing only
        the provided ``virt_obj`` nodes.
    virt_obj : list of str
        Sequence of virtual object identifiers to be added as nodes to the graph.
        Nodes that also appear in parsed rules are not duplicated.

    Returns
    -------
    networkx.DiGraph
        A directed graph with nodes for all operands, results and virtual objects.
        For each rule ``"A + B = C"`` there will be edges ``A -> C`` and ``B -> C``.

    Raises
    ------
    ValueError
        If a rule does not contain exactly two operands and one result when split
        using the expected separators (e.g. malformed string).
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


def repair(s: Union[str, list[str]]) -> Tuple[List[List[str]], Dict[str, List[str]]]:
    """
    Iteratively replace the most frequent adjacent symbol pairs in a string with new non-terminal symbols.

    Scans the input string for adjacent symbol pairs, identifies pairs that occur more than
    once, and replaces the most frequent pair with a freshly introduced non-terminal
    (e.g. ``A1``, ``A2``). Replacements are applied repeatedly until no pair occurs
    more than once. Returns the resulting symbol sequence and a mapping of introduced
    non-terminals to the pair they replaced.

    Parameters
    ----------
    s : str or list of str
        Input string of symbols (each character treated as a symbol). The function
        operates on single-character symbols and produces new multi-character
        non-terminal symbols of the form ``A{n}``.

    Returns
    -------
    symbols : list of str
        Sequence of symbols after all replacements. Contains original terminal
        symbols and introduced non-terminal symbols (strings like ``'A1'``).
    productions : dict
        Mapping from introduced non-terminal symbol to the two-symbol list it
        replaces, e.g. ``{'A1': ['a', 'b']}``.

    Raises
    ------
    TypeError
        If ``s`` is not a string.
    """
    if isinstance(s, str):
        symbols: List[List[str]] = [list(s)]
    else:
        if not all(isinstance(subs, str) for subs in s):
            raise TypeError("Input must be a string or a list of strings.")
        symbols: List[List[str]] = [list(subs) for subs in s]

    # Input safety check
    lowercase = set(string.ascii_lowercase)
    for symbol in symbols:
        assert lowercase.issuperset(symbol), "Input string must consist of lowercase ASCII characters only."

    symbols, rules = repair_pairs(symbols, lambda i: f'A{i + 1}')
    productions: Dict[str, List[str]] = {nt: list(pair) for nt, pair in rules}

    return symbols, productions


def convert_to_cnf(start_symbols: Union[str, List[str]],
                   productions: Dict[str, List[str]]) -> Tuple[str, Dict[str, List[str]]]:
    """
    Convert a context-free grammar (CFG) to Chomsky Normal Form (CNF).

    Transforms the given productions by:
    - mapping lowercase terminal symbols to fresh terminal non-terminals of the form ``T_<symbol>``;
    - introducing new non-terminals ``N<n>`` to break right-hand sides longer than two into binary rules;
    - ensuring a single start non-terminal (``S``) whose right-hand side has length one or two.

    Parameters
    ----------
    start_symbols : Union[str, List[str]]
        The start symbol or start sequence of symbols for the grammar. Can be a
        string (each character treated as a symbol) or a list of symbol strings.
    productions : dict
        Mapping from non-terminal symbols (keys) to lists of symbols (values). Each value is a list
        representing the right-hand side of a production. Terminals are assumed to be lowercase ASCII
        characters; other strings are treated as non-terminals.

    Returns
    -------
    start_nts : str
        The new start non-terminal (conventionally ``S``).
    cnf_productions : dict
        A dictionary mapping non-terminals to CNF-compliant right-hand sides. Right-hand sides are
        either a single terminal (e.g. ``['T_a']``) or two non-terminals (e.g. ``['N1', 'B']``).
        Terminal mappings for each original terminal are included as productions (``'T_a': ['a']``).

    Raises
    ------
    TypeError
        If ``start_symbol`` is not a string or if ``productions`` is not a mapping of non-terminals to lists.
    ValueError
        If a production contains an empty right-hand side.

    """
    cnf_productions: Dict[str, List[str]] = {}
    new_nt_counter: int = 1

    if isinstance(start_symbols, str):
        start_symbols = [start_symbols]

    # Map terminals to new non-terminals
    terminals = {s for exp in productions.values() for s in exp if s in string.ascii_lowercase}
    for start_symbol in start_symbols:
        terminals.update(s for s in start_symbol if s in string.ascii_lowercase)
    terminal_map: Dict[str, str] = {t: f'T_{t}' for t in terminals}
    cnf_productions.update({nt: [t] for t, nt in terminal_map.items()})

    def replace_terminals(symbols: List[str]) -> List[str]:
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
    for idx, word in enumerate(start_symbols):
        word = replace_terminals(list(word))
        start_nt = 'S_' + str(idx)
        # Fold left into binary rules, as for long productions above.
        if len(word) > 2:
            head = word[0]
            for symbol in word[1:-1]:
                new_nt = f'N{new_nt_counter}'
                cnf_productions[new_nt] = [head, symbol]
                head = new_nt
                new_nt_counter += 1
            word = [head, word[-1]]
        cnf_productions[start_nt] = word

    return start_nt, cnf_productions


def ai_core(s: Union[str, List[str]],
            debug: bool = False) -> Tuple[int, Dict[str, List[str]]]:
    """
    Convert an input string into Chomsky Normal Form (CNF) and compute a production count.

    This function performs a two-stage transformation: it first derives a set of
    productions and an augmented start symbol using the `repair` routine, then
    transforms those productions into CNF using `convert_to_cnf`. It also returns
    a scalar measure (`ai_count`) derived from the intermediate start symbol and
    the number of introduced productions.

    Parameters
    ----------
    s : Union[str, List[str]]
        Input string of symbols to be processed. Each character is treated as a
        terminal symbol for the initial `repair` pass.
    debug : bool, optional
        If True, print intermediate values (start symbol, productions, CNF
        productions and their lengths) to stdout for debugging. Default is False.

    Returns
    -------
    ai_count : int
        Integer count computed as `len(start_symbol) - 1 + len(productions)` where
        `start_symbol` and `productions` are the outputs of `repair(s)`. Intended
        as a simple complexity/path-length metric used by downstream logic.
    cnf_productions : dict
        Mapping of non-terminal symbols to CNF-compliant right-hand sides. Each
        value is a list of one terminal-nonterminal mapping (e.g. ``['T_a']``) or
        two non-terminals (e.g. ``['N1', 'B']``).

    Raises
    ------
    TypeError
        If ``s`` is not a string or a list of strings.
    """
    start_symbols, productions = repair(s)
    start_nts, cnf_productions = convert_to_cnf(start_symbols, productions)
    if debug:
        print(f"Start symbols: {start_symbols}", flush=True)
        print(f"Productions: {productions}", flush=True)
        print(f"Length of Productions: {len(productions)}", flush=True)
        print(f"CNF Productions: {cnf_productions}", flush=True)
        print(f"Length of CNF Productions: {len(cnf_productions)}", flush=True)
    # return len(cnf_productions) - len(set(s)), cnf_productions

    # Use temp to count number of terminal symbols (ai = # of non-terminal producing cnf production rules)
    temp = ""
    for obj in set(s):
        temp += obj

    return len(cnf_productions) - len(set(temp)), cnf_productions


def get_rules(s: Union[str, List[str]],
              production: Dict[str, List[str]],
              f_print: bool = False) -> List[str]:
    """
    Generate join rules from productions by performing a topological sort.

    Performs a Kahn-style topological traversal over the dependency graph
    defined by `production` and the initial symbols in `s`. Builds string
    mappings for intermediate non-terminals and emits join rules of the form
    ``"<left> + <right> = <result>"`` when a non-terminal expands to two symbols.

    Parameters
    ----------
    s : Union[str, List[str]]
        Input string or list of strings being processed.
    production : dict
        Mapping from non-terminal symbol to its right-hand side as a list of
        symbols. Right-hand sides are expected to be length 1 or 2. Keys that
        appear in `production` are treated as nodes that depend on the symbols in
        their value list.
    f_print : bool, optional
        If True, print processing diagnostics (start symbols, joins) to stdout.
        Default is False.

    Returns
    -------
    rules : list of str
        List of join rules in topologically valid order. Each rule is formatted
        as ``"A + B = C"`` where ``A`` and ``B`` are the left-hand constituents and
        ``C`` is the computed result string for the non-terminal.

    Raises
    ------
    TypeError
        If ``s`` is not a string or ``production`` is not a mapping-like object.
    ValueError
        If the dependency graph contains a cycle, the returned rules may be
        incomplete; callers should validate acyclicity prior to calling if full
        coverage is required.
    """
    in_degrees: Dict[str, int] = collections.defaultdict(int)
    adj: Dict[str, List[str]] = collections.defaultdict(list)
    tmap: Dict[str, str] = {}

    # Initialize in-degrees and adjacency list. Joint inputs need to be
    # flattened so their terminal symbols, rather than the complete strings,
    # become the roots of the dependency graph.
    initial_symbols = (
        s if isinstance(s, str) else (symbol for word in s for symbol in word)
    )
    for c in initial_symbols:
        in_degrees[c] = 0
    for course, prereqs in production.items():
        for req in prereqs:
            in_degrees[course] += 1
            adj[req].append(course)

    # Start queue with symbols having zero in-degrees
    start_q: collections.deque[str] = collections.deque(
        symbol for symbol, ins in in_degrees.items() if ins == 0
    )
    if f_print:
        print(f"Processing {s}", flush=True)
        print(f"Start symbols: {', '.join(start_q)}", flush=True)
        print("Joins:", flush=True)

    # Perform topological sort
    rules: List[str] = []
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


def extract_virtual_objects(rules: List[str]) -> List[str]:
    """
    Extract virtual objects from a list of join rules, excluding the final rule's result.

    Parameters
    ----------
    rules : list of str
        Sequence of join rules formatted as ``"A + B = C"``. An empty sequence
        results in an empty list.

    Returns
    -------
    virt_objs : list of str
        Sorted list of unique object identifiers found in the rules, excluding the
        right-hand side (result) of the last rule. Sorting is by string length
        (ascending).
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


def repair_with_pathways(s: Union[str, List[str]], f_print: bool = False, debug: bool = False) -> Tuple[
    int, List[str], nx.DiGraph]:
    """
    Compute pathway information from an input string: path length, virtual objects, and a rules graph.

    Performs grammar repair and CNF conversion via `ai_core`, derives join rules with `get_rules`,
    extracts virtual objects with `extract_virtual_objects`, and builds a directed graph of rules
    with `rules_to_graph`.

    Parameters
    ----------
    s : Union[str, List[str]]
        Input string or list of strings to be processed. Each character is treated as an initial terminal symbol.
    f_print : bool, optional
        If True, print join-rule processing diagnostics and the computed rules. Default is False.
    debug : bool, optional
        If True, enable debug printing in the underlying `ai_core` stage. Default is False.

    Returns
    -------
    ai_count : int
        Integer path-length metric returned by `ai_core` (heuristic count of introduced productions
        and start-symbol length).
    virt_obj : list of str
        Sorted list of virtual object identifiers extracted from the join rules (excluding the final
        rule result). Sorted by string length ascending.
    rules_graph : networkx.DiGraph
        Directed graph representing join relationships. For each join rule ``"A + B = C"`` there
        are edges ``A -> C`` and ``B -> C``; virtual objects are added as isolated nodes when present.

    Raises
    ------
    TypeError
        If ``s`` is not a string. Underlying routines may raise additional errors for malformed
        productions or rules.
    """
    # Get the production rules and path length
    path_len, productions = ai_core(s, debug=debug)
    # Get the rules
    rules = get_rules(s, productions, f_print=f_print)
    # Extract virtual objects
    virt_obj = extract_virtual_objects(rules)
    return path_len, virt_obj, rules_to_graph(rules, virt_obj)
