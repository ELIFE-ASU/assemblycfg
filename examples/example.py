import CFG as at

if __name__ == '__main__':
    # Example 1
    s_inpt = "abracadabra"
    ai, path = at.ai_upper_with_pathways(s_inpt)

    # Example 2
    s_inpt = "aaaaaaa"
    ai, path = at.ai_upper_with_pathways(s_inpt)

    # Example 3
    s_inpt = "ababcdcd"
    ai, path = at.ai_upper_with_pathways(s_inpt)
