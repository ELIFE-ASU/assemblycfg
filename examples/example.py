import CFG

if __name__ == '__main__':
    # Example 1
    s_inpt = "abracadabra"
    CFG.ai_with_pathways(s_inpt, f_print=True)
    print(flush=True)

    # Example 2
    s_inpt = "aaaaaaa"
    CFG.ai_with_pathways(s_inpt, f_print=True)
    print(flush=True)

    # Example 3
    s_inpt = "ababcdcd"
    CFG.ai_with_pathways(s_inpt, f_print=True)
    print(flush=True)
