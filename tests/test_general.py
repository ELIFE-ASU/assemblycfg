import CFG as at


def test_version():
    assert at.__version__ == "0.1.00"


def test_abracadabra():
    s_inpt = "abracadabra"
    ai, _ = at.ai_upper(s_inpt)
    ai_ref = 7
    assert ai == ai_ref


def test_path():
    s_inpts = ["abracadabra", "aaaaaaa", "ababcdcd"]
    ai_refs = [7, 4, 5]
    for i in range(len(s_inpts)):
        ai, path = at.ai_upper_with_pathways(s_inpts[i], f_print=False)

        final_path = path[-1].split(" = ")[1].strip()

        assert ai == ai_refs[i]
        assert len(path) == ai_refs[i]
        assert final_path == s_inpts[i]
