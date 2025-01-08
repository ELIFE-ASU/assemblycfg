import CFG
import matplotlib.pyplot as plt
import networkx as nx

if __name__ == '__main__':
    f_print = True

    # Example 1
    s_inpt = "abracadabra"
    ai, virt_obj, path = CFG.ai_with_pathways(s_inpt, f_print=f_print)
    print(f"AI: {ai}")
    print(f"virt_obj: {virt_obj}")
    # Plot the graph
    nx.draw(path, with_labels=True, font_weight='bold')
    plt.show()
    print(flush=True)

    # Example 2
    s_inpt = "aaaaaaa"
    CFG.ai_with_pathways(s_inpt, f_print=f_print)
    print(flush=True)

    # Example 3
    s_inpt = "ababcdcd"
    CFG.ai_with_pathways(s_inpt, f_print=f_print)
    print(flush=True)
