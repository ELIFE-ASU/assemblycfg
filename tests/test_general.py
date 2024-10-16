import CFG as at


def test_version():
    assert at.__version__ == '0.0.01'


def test_abracadabra():
    s_inpt = "abracadabra"
    ai, _ = at.ai_upper(s_inpt)
    ai_ref = 7
    assert ai == ai_ref
