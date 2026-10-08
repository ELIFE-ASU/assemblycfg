# Context-Free Grammars and String Assembly Index

[![PyPI](https://img.shields.io/pypi/v/assemblycfg.svg)](https://pypi.org/project/assemblycfg/)
[![Python](https://img.shields.io/pypi/pyversions/assemblycfg.svg)](https://pypi.org/project/assemblycfg/)
[![CI](https://github.com/ELIFE-ASU/assemblycfg/actions/workflows/ci.yml/badge.svg)](https://github.com/ELIFE-ASU/assemblycfg/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20562899.svg)](https://doi.org/10.5281/zenodo.20562899)

`assemblycfg` places bounds on the assembly index of directed strings and
molecules. Upper bounds come from the RePair smallest-grammar algorithm, which
quickly finds a short assembly pathway but does not guarantee the shortest.
Lower bounds come from LZ factorisation and vector addition chains.

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

## Upper bounds

### Strings

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

### Molecules

`calculate_assembly_path_graph_repair` runs RePair on the molecular graph
itself, adapted from the GraphRePair-inspired bound in
[parallelassemblycpp](https://github.com/ELIFE-ASU/parallelassemblycpp). It
repeatedly joins the most common pair of incident fragments, counting the most
occurrences that share no bond, so branched and cyclic motifs can be reused as
well as chains:

```python
import assemblycfg as cfg

smiles = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"
molgraph = cfg.smi_to_nx(smiles)
length, virtual_objects, path = cfg.calculate_assembly_path_graph_repair(molgraph)
print(f"a(Cholesterol) <= {length}")
```

The virtual objects are NetworkX graphs of the molecular fragments. Ties are
broken in atom order, so the default is deterministic; `iterations` adds passes
over random bond orders and keeps the shortest pathway. `graph_repair` returns
the full construction certificate. More complete programs are available in the
[`examples`](https://github.com/ELIFE-ASU/assemblycfg/tree/main/examples)
directory.

## Lower bounds

### LZ factorisation

`lz_lower_bound` places a valid lower bound on the assembly index of a single
directed string:

```python
import assemblycfg as cfg

print(cfg.lz_lower_bound("abracadabra"))  # 7
```

Read left to right, every assembly pathway builds the string by appending one
character or a substring that already occurs earlier in it. The fewest such
steps, found by dynamic programming over prefixes, bounds the assembly index
from below.

### Vector addition chains

`vac_lower_bound` places a valid lower bound on the assembly index of a string,
a molecule, or a list of them (bounding their joint assembly index):

```python
import assemblycfg as cfg

print(cfg.vac_lower_bound("abracadabra"))          # 7
print(cfg.vac_lower_bound(cfg.smi_to_nx("CCO")))   # 6
```

Counting the copies of each basic unit (characters, or bonds by element pair
and bond order) maps an object to a vector, and joining two objects adds their
vectors. Every assembly pathway is therefore a vector addition chain of the
same length, so the shortest chain reaching the object's vector bounds its
assembly index from below.

Chains are solved by [`vac`](https://github.com/ELIFE-ASU/additionchains),
which is installed on first use with `cargo` (install Rust from
[rustup.rs](https://rustup.rs)) into `~/.cache/assemblycfg/vac`. Set `VAC_PATH`
to use your own build, or `VAC_REF` to install a specific branch, tag or
commit. Without `vac`, the function warns and returns the closed-form bounds
it would start from.

## Package layout

| Module | Contents |
| --- | --- |
| `string_repair` | RePair upper bounds on string assembly index |
| `molecule_repair` | Graph RePair upper bounds on molecular assembly index |
| `lz` | LZ factorisation lower bound for strings |
| `vac` | Vector addition chain lower bounds for strings and molecules |
| `molecules` | Conversion between SMILES, Molfiles, RDKit and NetworkX |

Every public function is also available from the top-level `assemblycfg`
namespace.

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
