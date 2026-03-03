import matplotlib.pyplot as plt
import networkx as nx

import assemblycfg

if __name__ == '__main__':
    f_print = True

    # Example 1
    s_inpt = "abracadabra"
    path_len, virt_obj, path = assemblycfg.repair_with_pathways(s_inpt, f_print=f_print)
    print(f"assembly index =< {path_len}")
    print(f"virt_obj: {virt_obj}")
    # Plot the graph
    nx.draw(path, with_labels=True, font_weight='bold')
    plt.show()
    print(flush=True)

    # Example 2
    s_inpt = "aaaaaaa"
    path_len, virt_obj, path = assemblycfg.repair_with_pathways(s_inpt, f_print=f_print)
    print(f"assembly index =< {path_len}")
    print(f"virt_obj: {virt_obj}")
    print(flush=True)

    # Example 3
    s_inpt = "ababcdcd"
    path_len, virt_obj, path = assemblycfg.repair_with_pathways(s_inpt, f_print=f_print)
    print(f"assembly index =< {path_len}")
    print(f"virt_obj: {virt_obj}")
    print(flush=True)
