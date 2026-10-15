import csv
import os
from string import ascii_lowercase

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.container import ErrorbarContainer
from matplotlib.ticker import NullLocator, PercentFormatter
from rdkit import Chem

import assemblycfg as CFG
from random_strings import best_of
from random_strings_tradeoff import panel_grid, point

# Setting plot aesthetics for better visibility
plt.rcParams['axes.linewidth'] = 2.0

# The 20 standard amino acids, as one-letter codes
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
# Peptide lengths in residues, one panel each, in rows of PANELS_PER_ROW
LENGTHS = [4, 6, 8, 12, 20]
PANELS_PER_ROW = 3
N_SAMPLES = 10
SEED = 2024
# AssemblyCPP timeouts in seconds. A run stopped by the timeout returns the best upper bound found.
TIMEOUTS = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100]
# AssemblyCPP results at each timeout, written by an offline script
TABLE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "peptide_tradeoff.csv")


def random_peptides(length: int, n_samples: int = N_SAMPLES, seed: int = SEED) -> list:
    """
    Sample random peptides with residues drawn uniformly from the 20 standard amino acids.

    The peptides depend only on the length and seed, so the lookup table and this
    example draw the same peptides.

    Parameters:
        length (int): The number of residues in each peptide.
        n_samples (int): The number of peptides to sample.
        seed (int): The random seed.

    Returns:
        list of str: The SMILES strings of the sampled peptides.
    """
    rng = np.random.default_rng([seed, length])
    return [Chem.MolToSmiles(Chem.MolFromSequence("".join(seq)))
            for seq in rng.choice(AMINO_ACIDS, size=(n_samples, length))]


def heavy_atom_graph(smiles: str):
    """The heavy-atom molecular graph that every method is given."""
    return CFG.smi_to_nx(smiles, add_hydrogens=False)


def load_table(file: str = TABLE_FILE) -> dict:
    """
    Load the AssemblyCPP timeout lookup table.

    Runs at timeouts beyond one that finished exactly, or that ran out of memory, are
    not in the table, and take the result of that run.

    Parameters:
        file (str): Path to the CSV lookup table.

    Returns:
        dict: Maps each (length, sample) to a dict of timeout -> (assembly index, time, exact).
              The assembly index is an upper bound when inexact, and -1 when no bound was found.
    """
    table = {}
    with open(file, newline="") as f:
        for row in csv.DictReader(f):
            runs = table.setdefault((int(row["length"]), int(row["sample"])), {})
            runs[float(row["timeout_s"])] = (int(row["assembly_index"]), float(row["time_s"]),
                                             row["exact"] == "True", row["out_of_memory"] == "True")
    for key, runs in table.items():
        last = None
        for timeout in TIMEOUTS:
            if timeout in runs:
                last = runs[timeout]
            elif last is not None and (last[2] or last[3]):
                runs[timeout] = last
        table[key] = {timeout: run[:3] for timeout, run in runs.items()}
    return table


def graph_repair_upper_bound(graph) -> int:
    """Upper bound on the assembly index from graph RePair."""
    return CFG.graph_repair(graph, iterations=1, seed=0).upper_bound


def vac_closed_form_bound(graph) -> int:
    """Lower bound on the assembly index from vac's closed-form bounds, without its chain search."""
    return CFG.vac_lower_bound(graph, use_vac=False)


if __name__ == "__main__":
    table = load_table()

    n_rows = -(-len(LENGTHS) // PANELS_PER_ROW)
    fig = plt.figure(figsize=(4.4 * PANELS_PER_ROW, 4.3 * n_rows + 0.6))
    panels = panel_grid(fig, len(LENGTHS), PANELS_PER_ROW)
    xs = 16

    for (ax, row_start), length, letter in zip(panels, LENGTHS, ascii_lowercase):
        graphs = [heavy_atom_graph(smiles) for smiles in random_peptides(length)]
        runs = [table[length, sample] for sample in range(N_SAMPLES)]

        repair = [best_of(graph_repair_upper_bound, g) for g in graphs]
        vac = [best_of(vac_closed_form_bound, g) for g in graphs]

        # The exact assembly index where AssemblyCPP finished, else the best upper bound found
        n_exact = sum(any(e for _, _, e in r.values()) for r in runs)
        reference = np.array([min([ai for ai, _, _ in r.values() if ai >= 0] + [rp])
                              for r, (rp, _) in zip(runs, repair)])

        def error(bounds):
            return (np.asarray(bounds) - reference) / reference

        if n_exact < N_SAMPLES:
            # Where the true assembly index lies: at the reference when exact, else between the
            # VAC lower bound and the best upper bound
            lower = [ref if any(e for _, _, e in r.values()) else b
                     for r, ref, (b, _) in zip(runs, reference, vac)]
            ax.axhspan(error(lower).mean(), 0, color='grey', alpha=0.15, lw=0, label='Range of True Value')
        ax.axhline(0, color='grey', lw=1)

        # AssemblyCPP's best bound within each timeout. Its search is deterministic, so a bound found
        # within a shorter timeout is also found within a longer one. Before any bound is found the
        # trivial bound, one join per bond after the first, holds. Runs stopped by the timeout or by
        # running out of memory can overrun the timeout, so a run's time is capped at its timeout.
        trivial = [g.number_of_edges() - 1 for g in graphs]
        acpp_time = [np.mean([min(r[timeout][1], timeout) for r in runs]) for timeout in TIMEOUTS]
        acpp_err = []
        for i in range(len(TIMEOUTS)):
            best = [min([r[t][0] for t in TIMEOUTS[:i + 1] if r[t][0] >= 0] + [trivial_bound])
                    for r, trivial_bound in zip(runs, trivial)]
            acpp_err.append(error(best))
        acpp_err = np.array(acpp_err)
        ax.fill_between(acpp_time, acpp_err.min(axis=1), acpp_err.max(axis=1), color='green', alpha=0.25, lw=0)
        ax.plot(acpp_time, acpp_err.mean(axis=1), '-o', color='green', lw=2.5, ms=6, label='AssemblyCPP')

        point(ax, [t for _, t in repair], error([u for u, _ in repair]), 's', 'Graph RePair Upper Bound')
        point(ax, [t for _, t in vac], error([v for v, _ in vac]), 'v', 'Vector Chain Closed-Form Lower Bound')

        ax.set_xscale('log')
        # Linear within 10% of the reference, where the bounds are compared, and logarithmic beyond
        ax.set_yscale('symlog', linthresh=0.1, linscale=2)
        ax.set_yticks([-0.5, -0.05, 0, 0.05, 0.1, 1])
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.yaxis.set_minor_locator(NullLocator())
        bonds = np.mean([g.number_of_edges() for g in graphs])
        ax.set_title(f"{length} Residues, ~{bonds:.0f} Bonds ({n_exact}/{N_SAMPLES} Exact)", fontsize=xs - 2)
        ax.text(0.04, 0.96, f"({letter})", transform=ax.transAxes, ha='left', va='top',
                fontsize=xs, fontweight='bold')
        ax.set_xlabel('Compute Time (s)', fontsize=xs - 2)
        ax.tick_params(axis='both', which='major', labelsize=xs - 2, direction='in', length=6, width=2)
        ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1)
        ax.tick_params(axis='both', which='both', top=True, right=True)
        if not row_start:
            ax.tick_params(axis='y', labelleft=False)

    fig.supylabel('Relative Error in Joining Operations', fontsize=xs - 2)
    # The legend shows the markers of the bounds without their whiskers
    handles, labels = panels[-1][0].get_legend_handles_labels()
    handles = [h[0] if isinstance(h, ErrorbarContainer) else h for h in handles]
    fig.legend(handles, labels, loc='lower center', ncol=len(labels), fontsize=xs - 3, frameon=False)

    # Spacing set by hand, since tight_layout would widen every gap to fit the y labels of an
    # offset row's first panel
    fig.subplots_adjust(left=0.085, right=0.99, top=0.96, bottom=0.13, wspace=0.3, hspace=0.35)
    plt.savefig('peptide_tradeoff.png', dpi=600)
    plt.savefig('peptide_tradeoff.pdf')
    plt.show()
