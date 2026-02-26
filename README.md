# Context-Free Grammar (CFG) string assembly

Directed string assembly index calculator using the smallest grammar algorithm re-pair. This will find a short, valid path.
But there is no guarantee that it will find the shortest possible valid path.

CFG/cfg_ai.py has two useful functions: ai_upper and ai_upper_with_pathways. Both return the same path length, but
ai_upper_with_pathways prints the joining operations of the path.

## Installation

Prerequisites: 
networkx >= 3.4.2
rdkit >=2024.03.5

```
conda install conda-forge::networkx
```

Use pip to install this package.

```
pip install git+https://github.com/ELIFE-ASU/assemblycfg.git
```

## Examples

See the examples folder for examples of how to use the package.


