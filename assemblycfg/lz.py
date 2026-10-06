"""
Lower bound on string assembly index from an LZ-style factorisation.

An assembly pathway for a directed string can be replayed left to right: each
step appends either one character or a substring that already occurs earlier
in the string. The fewest such steps needed to build the string, found by
dynamic programming over its prefixes, is a lower bound on its assembly index.
"""

from typing import Dict, List, Tuple


def lz_lower_bound(string: str) -> int:
    """
    Bound the assembly index of a directed string from below.

    ``dp[i]``, the bound for ``string[:i + 1]``, extends a shorter prefix by
    one character, or by a suffix ``string[j:i + 1]`` that occurs in
    ``string[:j]``, at the cost of one step either way.

    Parameters
    ----------
    string : str
        The target string. Its first character is free.

    Returns
    -------
    int
        The lower bound.

    Raises
    ------
    TypeError
        If ``string`` is not a string.
    ValueError
        If ``string`` is empty.

    Notes
    -----
    ``dp`` never decreases, so the cheapest suffix to reuse is the longest.
    The start ``j`` of that suffix never decreases with ``i`` either, so one
    pass slides the window ``string[j:i + 1]`` across a suffix automaton
    that records where each substring first ends. The bound takes linear
    time.

    Examples
    --------
    >>> import assemblycfg as cfg
    >>> cfg.lz_lower_bound("abracadabra")
    7
    """
    if not isinstance(string, str):
        raise TypeError(f"Input must be a string, not {type(string).__name__}")
    if not string:
        raise ValueError("Input must be a non-empty string")

    end, link, length, edges = _suffix_automaton(string)
    dp = [0] * len(string)
    # The window string[j:i + 1] is the automaton state v, of size n.
    j, v, n = 1, 0, 0
    for i in range(1, len(string)):
        v = edges[v][string[i]]
        n += 1
        # Shrink the window until its first occurrence ends before j.
        while n and end[v] >= j:
            j += 1
            n -= 1
            if n == length[link[v]]:
                v = link[v]
        dp[i] = dp[j - 1] + 1 if n else dp[i - 1] + 1
    return dp[-1]


def _suffix_automaton(string: str) -> Tuple[List[int], List[int], List[int], List[Dict[str, int]]]:
    """
    Build the suffix automaton of *string*.

    Returns, per state, the end of the first occurrence of its substrings,
    its suffix link, the length of its longest substring and its edges.
    State 0 is the empty string.
    """
    end, link, length, edges = [-1], [-1], [0], [{}]
    last = 0
    for i, char in enumerate(string):
        cur = len(length)
        end.append(i)
        link.append(0)
        length.append(length[last] + 1)
        edges.append({})
        p = last
        while p != -1 and char not in edges[p]:
            edges[p][char] = cur
            p = link[p]
        if p != -1:
            q = edges[p][char]
            if length[p] + 1 == length[q]:
                link[cur] = q
            else:
                clone = len(length)
                end.append(end[q])
                link.append(link[q])
                length.append(length[p] + 1)
                edges.append(dict(edges[q]))
                while p != -1 and edges[p].get(char) == q:
                    edges[p][char] = clone
                    p = link[p]
                link[q] = link[cur] = clone
        last = cur
    return end, link, length, edges
