"""
Lower bound on string assembly index from an LZ-style factorisation.

An assembly pathway for a directed string can be replayed left to right: each
step appends either one character or a substring that already occurs earlier
in the string. The fewest such steps needed to build the string, found by
dynamic programming over its prefixes, is a lower bound on its assembly index.
"""


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
    ``dp`` never decreases, so the cheapest suffix to reuse is the longest,
    and whether ``string[j:i + 1]`` occurs in ``string[:j]`` is monotone in
    ``j``. A binary search over ``j`` therefore finds it.

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

    dp = [0] * len(string)
    for i in range(1, len(string)):
        dp[i] = dp[i - 1] + 1
        # Smallest j > 0 whose suffix string[j:i + 1] occurs in string[:j].
        low, high = 1, i
        while low <= high:
            mid = (low + high) // 2
            if string[mid:i + 1] in string[:mid]:
                dp[i] = min(dp[i], dp[mid - 1] + 1)
                high = mid - 1
            else:
                low = mid + 1
    return dp[-1]
