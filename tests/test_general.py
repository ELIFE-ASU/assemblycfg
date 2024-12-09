import CFG


def test_abracadabra():
    # Test input string
    s_inpt = "abracadabra"
    # Get the assembly calculator result
    ai, _ = CFG.ai_upper(s_inpt)
    # Assembly calculator result gives 7
    ai_ref = 7
    # Check if the result is correct
    assert ai == ai_ref


def test_path():
    # Test input strings
    s_inpts = ["abracadabra", "aaaaaaa", "ababcdcd"]
    # Assembly calculator results
    ai_refs = [7, 4, 5]
    for i in range(len(s_inpts)):
        # Get the assembly calculator result
        ai, path = CFG.ai_upper_with_pathways(s_inpts[i], f_print=False)
        # Get the final path
        final_path = path[-1].split(" = ")[1].strip()
        # Check if the result is correct
        assert ai == ai_refs[i]
        assert len(path) == ai_refs[i]
        assert final_path == s_inpts[i]
