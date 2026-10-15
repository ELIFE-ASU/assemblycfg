# Examples

Run each script from its own folder; figures are written next to it. The
plotting scripts need the `plot` extra (`pip install "assemblycfg[plot]"`).

## Usage

| Script | Shows |
| --- | --- |
| `strings/pathways.py` | RePair upper bounds and assembly pathways for short strings |
| `molecules/lipids.py` | Graph RePair upper bounds and vector addition chain lower bounds on fatty acids, triglycerides and phospholipids |
| `molecules/gdb_bounds.py` | Graph RePair upper bounds and vector addition chain lower bounds against the exact assembly index on one million GDB-17 molecules |

## Benchmarks

Accuracy and resource use of the string bounds on random strings over a
four-letter alphabet. The exact assembly indices come from AssemblyCPP, stored
in the CSV tables by offline scripts.

| Script | Shows |
| --- | --- |
| `benchmarks/random_strings.py` | RePair and LZ bounds against the exact assembly index, and compute time against string length |
| `benchmarks/random_strings_tradeoff.py` | Error against compute time for RePair, LZ and AssemblyCPP at increasing timeouts |
| `benchmarks/random_strings_memory_scaling.py` | Peak memory of RePair and LZ against string length |

The molecular benchmarks use graph RePair and the closed-form vector addition
chain bound, without vac's chain search. The tradeoff benchmarks compare them
with AssemblyCPP on random peptides, and on GDB-17 molecules from
`molecules/1M-GDB-Combined.csv`, whose exact assembly indices are the reference.

| Script | Shows |
| --- | --- |
| `benchmarks/peptide_tradeoff.py` | Error against compute time for graph RePair, the closed-form VAC bound and AssemblyCPP at increasing timeouts, by peptide length |
| `benchmarks/gdb_tradeoff.py` | The same comparison on GDB-17 molecules, by bond count |
| `benchmarks/peptide_memory_scaling.py` | Peak memory of graph RePair and the closed-form VAC bound against bond count, on random peptides |
