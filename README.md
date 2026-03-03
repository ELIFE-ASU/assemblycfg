# Context-Free Grammars and String Assembly Index

Directed string assembly index calculator using the smallest grammar algorithm RePair. This will quickly find a short assembly path, but there is no guarantee that it will find the shortest possible assembly path. Thus, this path length serves as an upper bound to the assembly index. This method works best on strings but can also be applied to molecular graphs as we will demonstrate below.

## Installation

Prerequisites: 
networkx >= 3.4.2
rdkit >=2024.03.5
matplotlib>=3.9.2

Use pip to install this package.

```
pip install assemblycfg
```

## Examples

The central function of this package, `cfg.repair_with_pathways` returns three items. First it returns the integer path length with upper bounds the assembly index, second it returns the list of virtual object strings which were used along the assembly path identified by RePair, and third it returns a networkx DiGraph object depicting the assembly path.

```
import assemblycfg as cfg
l, vo, path = cfg.repair_with_pathways("abracadabra")
print(f'a("abracadabra") =< {l}')
print(f"Virtual objects used: {vo}")
```
You can visualize the pathway as follows
```
import networkx as nx
import matplotlib.pyplot as plt
nx.draw(path, with_labels=True, font_weight='bold', pos=nx.spring_layout(path))
plt.show()
```
though these pathway visuals easy get unweildy. We recommend the python package AssemblyTheoryTools for more sophisticated pathway plotting functions.

One can also apply these methods to molecular assembly index. The function `calculate_assembly_path_det` can place a valid upper bound on the assembly index of any molecule, though it performs best on 'stringy' molecules like lipids. Starting from a SMILES string for cholesterol, we convert it into a networkx graph format before passing it to the calculator.
```
import assemblycfg as cfg
smi_str = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C" # SMILES string for cholesterol
molgraph = cfg.smi_to_nx(smi_to_nx)
l, vo, path = cfg.calculate_assembly_path_det(molgraph)
print(f'a(Cholesterol) =< {l}')
```
These virtual objects will also be networkx graphs representing molecular fragments.

See the examples folder for more examples of how to use the package.

These algorithms are described in Siebert et al. (In Prep); if you find this package useful, please cite this paper.
