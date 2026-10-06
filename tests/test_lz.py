import pytest

import assemblycfg as cfg


@pytest.mark.parametrize(
    ("string", "bound"),
    [
        ("a", 0),
        ("ab", 1),
        ("aaaaaaaa", 3),
        ("abab", 2),
        ("abcd", 3),
        ("abracadabra", 7),
    ],
)
def test_lz_lower_bound_known_values(string, bound):
    assert cfg.lz_lower_bound(string) == bound


@pytest.mark.parametrize("string", ["abracadabra", "abcabcabcabc", "aabbaabbab"])
def test_lz_lower_bound_does_not_exceed_repair_upper_bound(string):
    upper, _, _ = cfg.repair_with_pathways(string)
    assert cfg.lz_lower_bound(string) <= upper


def test_lz_lower_bound_rejects_empty_string():
    with pytest.raises(ValueError):
        cfg.lz_lower_bound("")


def test_lz_lower_bound_rejects_non_strings():
    with pytest.raises(TypeError):
        cfg.lz_lower_bound(["ab", "cd"])
