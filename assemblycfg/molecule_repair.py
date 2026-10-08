"""
Upper bounds on molecular assembly index from RePair on the graph itself.

Ported from ``graphRepair.h`` in parallelassemblycpp. Active tokens partition
the bonds into connected fragments. A production joins two incident,
edge-disjoint tokens. Exact isomorphism classes of the full expanded labelled
fragment identify reusable products, including joins at several vertices. All
vertices remain available as attachment sites; this deliberately relaxes
GraphRePair's fixed boundary interfaces to match molecular assembly. It is not
the paper's linear-time compressor.

A new binary production costs one join and reusing an existing one costs zero,
so replacing k disjoint pairs saves k - 1 (or k for an existing product)
against the residual join count. At termination, r connected tokens covering
c nonempty components can be joined in r - c steps. The default disconnected
convention adds c - 1, giving ``rules + r - 1``; ``compensate_disjoint`` omits
those c - 1 steps. Isolated atoms are not bond primitives.

Two changes reduce the C++ version's dependence on bond order. A pair's count
is the maximum number of its occurrences that share no token, found by
maximum matching rather than greedily in bond order. And further passes may
retry the compression on random depth-first bond orders, keeping the best.
"""
import random
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
from networkx.algorithms.isomorphism import categorical_edge_match, categorical_node_match

_NODE_MATCH = categorical_node_match('color', None)
_EDGE_MATCH = categorical_edge_match('color', None)


@dataclass
class Fragment:
    """A connected fragment: its symbol and sorted indices into ``edges``."""
    symbol: int
    edges: Tuple[int, ...]


@dataclass
class Rule:
    """Binary production ``id -> left right`` with one witnessing occurrence."""
    id: int
    left: int
    right: int
    left_edges: Tuple[int, ...]
    right_edges: Tuple[int, ...]
    edges: Tuple[int, ...]


@dataclass
class GraphRepairResult:
    """A replayable certificate: fragment edges index into ``edges``."""
    upper_bound: int = 0
    trivial_upper_bound: int = 0
    rule_count: int = 0
    remaining_fragments: int = 0
    components: int = 0
    compensate_disjoint: bool = False
    edges: List[Tuple[Any, Any]] = field(default_factory=list)
    terminals: List[Fragment] = field(default_factory=list)
    rules: List[Rule] = field(default_factory=list)
    residual: List[Fragment] = field(default_factory=list)


class _IsomorphismClasses:
    """
    Intern labelled graphs by exact isomorphism class.

    Graphs are bucketed by a Weisfeiler-Lehman hash and confirmed against each
    bucket's representatives with VF2, so distinct classes never share an ID.
    """

    def __init__(self):
        self._buckets: Dict[Tuple[int, int, str], List[Tuple[nx.Graph, int]]] = {}
        self._count = 0

    def __call__(self, graph: nx.Graph) -> int:
        key = (graph.number_of_nodes(), graph.number_of_edges(),
               nx.weisfeiler_lehman_graph_hash(graph, node_attr='label', edge_attr='label'))
        bucket = self._buckets.setdefault(key, [])
        for representative, class_id in bucket:
            if nx.is_isomorphic(graph, representative, node_match=_NODE_MATCH, edge_match=_EDGE_MATCH):
                return class_id
        class_id = self._count
        self._count += 1
        bucket.append((graph, class_id))
        return class_id


def _edge_list(graph: nx.Graph) -> List[Tuple[Any, Any, Any]]:
    """Edges in node order, each from its earlier endpoint, as in ``molGraph::writeEdgeList``."""
    position = {node: i for i, node in enumerate(graph.nodes)}
    return [(u, v, data['color'])
            for u in graph.nodes
            for v, data in graph.adj[u].items()
            if position[u] < position[v]]


def _fragment(graph: nx.Graph, edges: List[Tuple[Any, Any, Any]], selected) -> nx.Graph:
    """The labelled subgraph spanned by the selected edges, with WL hash labels."""
    fragment = nx.Graph()
    for i in selected:
        u, v, color = edges[i]
        for node in (u, v):
            if node not in fragment:
                color_u = graph.nodes[node]['color']
                fragment.add_node(node, color=color_u, label=repr(color_u))
        fragment.add_edge(u, v, color=color, label=repr(color))
    return fragment


def _merge(left: Tuple[int, ...], right: Tuple[int, ...]) -> Tuple[int, ...]:
    return tuple(sorted(set(left) | set(right)))


def _random_edge_list(graph: nx.Graph, rng: random.Random) -> List[Tuple[Any, Any, Any]]:
    """Edges in a random depth-first order, so bonded fragments stay close in the list."""
    shuffled = nx.Graph()
    shuffled.add_nodes_from(rng.sample(list(graph), len(graph)))
    for u in shuffled:
        for v in rng.sample(list(graph.adj[u]), len(graph.adj[u])):
            shuffled.add_edge(u, v)
    return [(u, v, graph.edges[u, v]['color']) for u, v in nx.edge_dfs(shuffled, list(shuffled))]


def _disjoint(occurrences: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    The most occurrences that share no token, in their original order.

    Occurrences may share vertices but not tokens (and hence bonds), so this
    is a maximum matching in the graph whose edges are the occurrences.
    """
    used = set()
    selected = []
    for first, second in occurrences:
        if first in used or second in used:
            continue
        used.update((first, second))
        selected.append((first, second))
    tokens = {token for occurrence in occurrences for token in occurrence}
    if len(selected) == len(tokens) // 2 or len(selected) == len(occurrences):
        return selected  # greedy is already maximum
    conflicts = nx.Graph(occurrences)
    matching = {frozenset(pair) for pair in nx.max_weight_matching(conflicts, maxcardinality=True)}
    return [occurrence for occurrence in occurrences if frozenset(occurrence) in matching]


def _compress(graph: nx.Graph, edges: List[Tuple[Any, Any, Any]], compensate_disjoint: bool) -> GraphRepairResult:
    """One compression pass; ``edges`` order breaks every tie."""
    result = GraphRepairResult(compensate_disjoint=compensate_disjoint,
                               edges=[(u, v) for u, v, _ in edges])
    bonded = graph.subgraph(node for node in graph if graph.degree(node) > 0)
    result.components = nx.number_connected_components(bonded)
    if not edges:
        return result
    adjustment = result.components if compensate_disjoint else 1
    result.trivial_upper_bound = len(edges) - adjustment

    classes = _IsomorphismClasses()
    dictionary: Dict[int, int] = {}  # isomorphism class -> symbol
    tokens: List[Tuple[Fragment, frozenset]] = []
    for i, (u, v, _) in enumerate(edges):
        class_id = classes(_fragment(graph, edges, (i,)))
        if class_id not in dictionary:
            dictionary[class_id] = len(dictionary)
            result.terminals.append(Fragment(dictionary[class_id], (i,)))
        tokens.append((Fragment(dictionary[class_id], (i,)), frozenset((u, v))))
    active = list(range(len(edges)))

    # Token IDs are stable, so unchanged adjacent pairs keep their classes.
    pair_classes: Dict[Tuple[int, int], int] = {}
    while len(active) > 1:
        rank = {token: i for i, token in enumerate(active)}
        at_vertex: Dict[Any, List[int]] = {}
        for token in active:
            for vertex in tokens[token][1]:
                at_vertex.setdefault(vertex, []).append(token)

        groups: Dict[int, List[Tuple[int, int]]] = {}  # class -> occurrences, first seen first
        for first in active:
            neighbours = {other for vertex in tokens[first][1] for other in at_vertex[vertex]
                          if rank[other] > rank[first]}
            for second in sorted(neighbours, key=rank.__getitem__):
                pair = (min(first, second), max(first, second))
                if pair not in pair_classes:
                    pair_classes[pair] = classes(_fragment(
                        graph, edges, _merge(tokens[first][0].edges, tokens[second][0].edges)))
                groups.setdefault(pair_classes[pair], []).append((first, second))

        best_saving, best_class, best_pairs = 0, None, []
        for class_id, occurrences in groups.items():
            selected = _disjoint(occurrences)
            saving = len(selected) - (0 if class_id in dictionary else 1)
            if saving > best_saving:
                best_saving, best_class, best_pairs = saving, class_id, selected
        if best_class is None:
            break

        if best_class not in dictionary:
            dictionary[best_class] = len(dictionary)
            left, right = tokens[best_pairs[0][0]][0], tokens[best_pairs[0][1]][0]
            result.rules.append(Rule(dictionary[best_class], left.symbol, right.symbol,
                                     left.edges, right.edges, _merge(left.edges, right.edges)))
        symbol = dictionary[best_class]

        removed = set()
        replacements = []
        for first, second in best_pairs:
            (left, left_vertices), (right, right_vertices) = tokens[first], tokens[second]
            removed.update((first, second))
            replacements.append(len(tokens))
            tokens.append((Fragment(symbol, _merge(left.edges, right.edges)), left_vertices | right_vertices))
        active = replacements + [token for token in active if token not in removed]
        active.sort(key=lambda token: tokens[token][0].edges[0])

    result.residual = [tokens[token][0] for token in active]
    result.rule_count = len(result.rules)
    result.remaining_fragments = len(result.residual)
    result.upper_bound = result.rule_count + result.remaining_fragments - adjustment
    return result


def graph_repair(graph: nx.Graph,
                 compensate_disjoint: bool = False,
                 iterations: int = 1,
                 seed: Optional[int] = None) -> GraphRepairResult:
    """
    Calculate a constructive assembly upper bound by graph pair replacement.

    Starting from one token per bond, each round groups incident token pairs
    by the isomorphism class of their union and replaces the most occurrences
    of the class that saves the most joins, until no replacement saves one.
    Each round revisits every incident pair but canonicalises only the new
    ones, so a pass costs time quadratic in the number of bonds plus one
    isomorphism test per new pair.

    Parameters
    ----------
    graph : networkx.Graph
        Undirected graph whose nodes and edges carry a ``'color'`` attribute.
    compensate_disjoint : bool, optional
        If True, do not charge the joins between components that contain
        bonds. Default is False.
    iterations : int, optional
        Number of passes. The first breaks ties in the graph's own node and
        adjacency order, which is deterministic; the rest use random
        depth-first bond orders. The shortest pathway is kept. Default is 1.
    seed : int, optional
        Seed for the random bond orders. Default is None.

    Returns
    -------
    GraphRepairResult
        The bound, the trivial bound (bonds minus one), and a certificate of
        terminals, rules and residual fragments.

    Raises
    ------
    ValueError
        If ``iterations`` is less than one.
    """
    if iterations < 1:
        raise ValueError("iterations must be at least 1")
    best = _compress(graph, _edge_list(graph), compensate_disjoint)
    rng = random.Random(seed)
    for _ in range(iterations - 1):
        result = _compress(graph, _random_edge_list(graph, rng), compensate_disjoint)
        if result.upper_bound < best.upper_bound:
            best = result
    return best


def calculate_assembly_path_graph_repair(graph: nx.Graph,
                                         iterations: int = 1,
                                         compensate_disjoint: bool = False,
                                         seed: Optional[int] = None) -> Tuple[
    int, Optional[List[nx.Graph]], Optional[nx.DiGraph]]:
    """
    Upper-bound the assembly index of a molecular graph with graph RePair.

    Fragments are general connected subgraphs, so repeated branched motifs
    and ring closures can be reused as well as chains. See ``graph_repair``.

    Parameters
    ----------
    graph : networkx.Graph
        Undirected graph whose nodes and edges carry a ``'color'`` attribute.
    iterations : int, optional
        Number of passes; all but the first use random bond orders. Default is 1.
    compensate_disjoint : bool, optional
        If True, do not charge the joins between components that contain
        bonds. Default is False.
    seed : int, optional
        Seed for the random bond orders. Default is None.

    Returns
    -------
    upper_bound : int
        The number of joins in the constructed assembly pathway.
    virtual_objects : list of networkx.Graph or None
        One subgraph of ``graph`` per discovered symbol (terminals, then
        rules), indexed by symbol. ``None`` for graphs without bonds.
    path : networkx.DiGraph or None
        Edges point from each rule's symbol to its two component symbols.
        ``None`` for graphs without bonds.
    """
    result = graph_repair(graph, compensate_disjoint=compensate_disjoint, iterations=iterations, seed=seed)
    if not result.edges:
        return 0, None, None
    path = nx.DiGraph()
    virtual_objects = []
    for fragment in result.terminals + result.rules:
        path.add_node(len(virtual_objects))
        virtual_objects.append(graph.edge_subgraph(result.edges[i] for i in fragment.edges).copy())
    for rule in result.rules:
        path.add_edge(rule.id, rule.left)
        path.add_edge(rule.id, rule.right)
    return result.upper_bound, virtual_objects, path
