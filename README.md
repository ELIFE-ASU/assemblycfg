# Context-Free Grammars and String Assembly Index

[![PyPI](https://img.shields.io/pypi/v/assemblycfg.svg)](https://pypi.org/project/assemblycfg/)
[![Python](https://img.shields.io/pypi/pyversions/assemblycfg.svg)](https://pypi.org/project/assemblycfg/)
[![CI](https://github.com/ELIFE-ASU/assemblycfg/actions/workflows/ci.yml/badge.svg)](https://github.com/ELIFE-ASU/assemblycfg/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20562899.svg)](https://doi.org/10.5281/zenodo.20562899)

`assemblycfg` calculates upper bounds on directed string and molecular assembly
indices using the RePair smallest-grammar algorithm. It quickly finds a short
assembly path, but it does not guarantee the shortest possible path.

## Installation

`assemblycfg` supports Python 3.12 and later. Install the package from PyPI; its
runtime dependencies are installed automatically:

```bash
python -m pip install assemblycfg
```

Plotting is optional. To run the visual examples, install the `plot` extra:

```bash
python -m pip install "assemblycfg[plot]"
```

## String assembly

The central function, `repair_with_pathways`, returns an upper bound on the
assembly index, the virtual objects used along the path, and a NetworkX directed
graph representing that path:

```python
import assemblycfg as cfg

length, virtual_objects, path = cfg.repair_with_pathways("abracadabra")
print(f'a("abracadabra") <= {length}')
print(f"Virtual objects used: {virtual_objects}")
```

Inputs may be a lowercase ASCII string or a list of such strings for a joint
assembly path.

With the optional plotting dependency installed, the path can be visualized as
follows:

```python
import matplotlib.pyplot as plt
import networkx as nx

nx.draw(path, with_labels=True, font_weight="bold", pos=nx.spring_layout(path))
plt.show()
```

These graphs can become unwieldy. The
[AssemblyTheoryTools](https://pypi.org/project/assemblytheorytools/) package
provides more sophisticated pathway plotting functions.

## Molecular assembly

`calculate_assembly_path_det` places a valid upper bound on the assembly index
of a molecule. It performs especially well on string-like molecules such as
lipids. For example:

```python
import assemblycfg as cfg

smiles = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"
molgraph = cfg.smi_to_nx(smiles)
length, virtual_objects, path = cfg.calculate_assembly_path_det(molgraph)
print(f"a(Cholesterol) <= {length}")
```

The virtual objects returned by the molecular workflow are NetworkX graphs that
represent molecular fragments. More complete programs are available in the
[`examples`](https://github.com/ELIFE-ASU/assemblycfg/tree/main/examples)
directory.

## Development

Development dependencies use the standardized dependency-groups table. With a
recent version of pip:

```bash
python -m pip install --upgrade pip
python -m pip install --group dev -e .
python -m pytest
python -m build
python scripts/check_dist.py
```

Maintainers should follow the version and Trusted Publishing checklist in
[`RELEASING.md`](https://github.com/ELIFE-ASU/assemblycfg/blob/main/RELEASING.md).

## Citation

If you use this package, cite the archived software release at
[doi:10.5281/zenodo.20562899](https://doi.org/10.5281/zenodo.20562899). Complete
software citation metadata is provided in
[`CITATION.cff`](https://github.com/ELIFE-ASU/assemblycfg/blob/main/CITATION.cff).
