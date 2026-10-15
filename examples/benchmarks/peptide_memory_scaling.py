from functools import partial

import matplotlib.pyplot as plt
import numpy as np

from peptide_tradeoff import graph_repair_upper_bound, heavy_atom_graph, random_peptides, vac_closed_form_bound
from random_strings_memory_scaling import band, exponent, peak_memory

# Setting plot aesthetics for better visibility
plt.rcParams['axes.linewidth'] = 2.0

# Peptide lengths in residues
LENGTHS = np.unique(np.geomspace(2, 200, 13).astype(int))
N_SAMPLES = 10

# Label, function and line style
METHODS = [
    ('Graph RePair Upper Bound', graph_repair_upper_bound, '--'),
    ('Vector Chain Closed-Form Lower Bound', vac_closed_form_bound, ':'),
]


if __name__ == "__main__":
    graphs = {length: [heavy_atom_graph(smiles) for smiles in random_peptides(length, N_SAMPLES)]
              for length in LENGTHS}
    bonds = np.array([[g.number_of_edges() for g in graphs[length]] for length in LENGTHS])
    mean_bonds = bonds.mean(axis=1)

    # Load vac's chain table before measuring, so the first measurement does not count it
    vac_closed_form_bound(graphs[LENGTHS[0]][0])

    results = []
    for label, func, style in METHODS:
        peaks = np.array([[peak_memory(func, g) for g in graphs[length]] for length in LENGTHS])
        k = exponent(mean_bonds, peaks)
        print(f"{label}: {peaks[-1].mean() / 1e6:.2f} MB at ~{mean_bonds[-1]:.0f} bonds, "
              f"peak memory ~ bonds^{k:.2f}", flush=True)
        results.append((label, k, style, peaks))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    xs = 16

    # The exponents describe the total peak memory, so only its legend shows them
    for label, k, style, peaks in results:
        band(ax1, mean_bonds, peaks / 1e6, style, rf'{label} ($\propto n^{{{k:.2f}}}$)')
        band(ax2, mean_bonds, peaks / bonds, style, label)

    ax1.set_title("Peak Memory vs Bond Count", fontsize=xs)
    ax1.set_ylabel('Peak Memory (MB)', fontsize=xs - 2)
    ax2.set_title("Peak Memory per Bond vs Bond Count", fontsize=xs)
    ax2.set_ylabel('Peak Memory per Bond (bytes)', fontsize=xs - 2)

    for ax in [ax1, ax2]:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(mean_bonds[0], mean_bonds[-1])
        ax.set_xlabel('Heavy-Atom Bond Count', fontsize=xs - 2)
        ax.legend(loc="upper left", fontsize=xs - 3)
        ax.tick_params(axis='both', which='major', labelsize=xs - 2, direction='in', length=6, width=2)
        ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1)
        ax.tick_params(axis='both', which='both', top=True, right=True)

    # Headroom above the curves for the legend, and on the per-bond axis, room
    # below down to a whole decade so that at least two decades are labelled
    ax1.set_ylim(top=ax1.get_ylim()[1] * 4)
    bottom, top = ax2.get_ylim()
    ax2.set_ylim(10 ** np.floor(np.log10(bottom)), top * 2)

    plt.tight_layout()
    plt.savefig('peptide_memory_scaling.png', dpi=600)
    plt.savefig('peptide_memory_scaling.pdf')
    plt.show()
