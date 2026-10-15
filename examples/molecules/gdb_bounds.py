import csv
import os
from multiprocessing import Pool

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import MaxNLocator

import assemblycfg as CFG

# Setting plot aesthetics for better visibility
plt.rcParams['axes.linewidth'] = 2.0

HERE = os.path.dirname(os.path.abspath(__file__))
# One million molecules sampled from GDB-17 with exact heavy-atom assembly indices from assemblytheorytools
TABLE_FILE = os.path.join(HERE, "1M-GDB-Combined.csv")
# Graph RePair upper bounds, cached because the full run takes a few CPU-hours
UPPER_CACHE_FILE = os.path.join(HERE, "gdb_graph_repair.npy")
# Vector addition chain lower bounds and whether vac proved each chain minimal
LOWER_CACHE_FILE = os.path.join(HERE, "gdb_vac_lower_bound.npy")
# Graph RePair passes per molecule; all but the first use random bond orders
ITERATIONS = 1
# Budget in seconds for each vac search; one that runs out still gives a valid, weaker bound
TIMEOUT = 2.0


def load_table(file: str = TABLE_FILE) -> tuple:
    """
    Load the SMILES strings and exact assembly indices.

    Parameters:
        file (str): Path to the CSV table with SMILES and Ai columns.

    Returns:
        tuple: A list of SMILES strings and an array of their assembly indices.
    """
    smiles, ai = [], []
    with open(file, newline="") as f:
        for row in csv.DictReader(f):
            smiles.append(row["SMILES"])
            ai.append(int(row["Ai"]))
    return smiles, np.array(ai)


def progress(results, total: int, label: str, every: int):
    """Yield from results, printing a count every so often."""
    for i, result in enumerate(results):
        if (i + 1) % every == 0 or i + 1 == total:
            print(f"{label}: {i + 1}/{total}", flush=True)
        yield result


def graph_repair_upper_bound(smiles: str) -> int:
    """Upper bound on the assembly index of the heavy-atom graph from graph RePair."""
    graph = CFG.smi_to_nx(smiles, add_hydrogens=False)
    return CFG.graph_repair(graph, iterations=ITERATIONS, seed=0).upper_bound


def upper_bounds(smiles: list, file: str = UPPER_CACHE_FILE) -> np.ndarray:
    """
    Calculate graph RePair upper bounds on every molecule across all cores.

    The bounds are saved to a cache file and loaded from it on later runs. Delete
    the file to recalculate them.

    Parameters:
        smiles (list of str): The SMILES strings.
        file (str): Path to the cache file.

    Returns:
        numpy.ndarray: The upper bound on each molecule's assembly index.
    """
    if os.path.exists(file):
        bounds = np.load(file)
        if len(bounds) == len(smiles):
            return bounds
    with Pool() as pool:
        bounds = np.array(list(progress(pool.imap(graph_repair_upper_bound, smiles, chunksize=1000),
                                        len(smiles), "graph RePair", 50000)), dtype=int)
    np.save(file, bounds)
    return bounds


def bond_counts(smiles: str) -> tuple:
    """
    Key a molecule by its bond-type counts, largest first.

    The vector addition chain bound depends only on these counts, so molecules
    with the same key share a bound. Disconnected molecules are keyed by their
    SMILES, since each component is a separate target.
    """
    graph = CFG.smi_to_nx(smiles, add_hydrogens=False)
    if not nx.is_connected(graph):
        return (smiles,)
    return tuple(sorted(CFG.mol_unit_counts(graph).values(), reverse=True))


def vac_bound(smiles: str) -> tuple:
    """Lower bound on the assembly index of the heavy-atom graph, and whether it is proven minimal."""
    bound, info = CFG.vac_lower_bound(CFG.smi_to_nx(smiles, add_hydrogens=False), timeout=TIMEOUT, return_info=True)
    return bound, info["proven_optimal"]


def lower_bounds(smiles: list, file: str = LOWER_CACHE_FILE) -> np.ndarray:
    """
    Calculate vector addition chain lower bounds on every molecule across all cores.

    vac runs once per distinct bond-count vector, which is far fewer than the
    molecules. The bounds are saved to a cache file and loaded from it on later
    runs. Delete the file to recalculate them.

    Parameters:
        smiles (list of str): The SMILES strings.
        file (str): Path to the cache file.

    Returns:
        numpy.ndarray: One row per molecule of the lower bound and whether vac proved it minimal.
    """
    if os.path.exists(file):
        bounds = np.load(file)
        if len(bounds) == len(smiles):
            return bounds
    with Pool() as pool:
        keys = list(progress(pool.imap(bond_counts, smiles, chunksize=1000), len(smiles), "bond counts", 50000))
        # One representative molecule per distinct vector
        first = {}
        for i, key in enumerate(keys):
            first.setdefault(key, i)
        unique = list(first)
        solved = list(progress(pool.imap(vac_bound, [smiles[first[key]] for key in unique], chunksize=1),
                               len(unique), "vac vectors", 1000))
    by_key = dict(zip(unique, solved))
    bounds = np.array([by_key[key] for key in keys], dtype=int)
    np.save(file, bounds)
    return bounds


def bound_panel(ax, ai, bound, norm, cmap, title, ylabel, xs):
    """
    Plot one point per (exact, bound) pair coloured by the number of molecules on it.

    Returns:
        matplotlib.collections.PathCollection: The scatter, for the colour bar.
    """
    pairs, counts = np.unique(np.column_stack([ai, bound]), axis=0, return_counts=True)
    order = np.argsort(counts)  # Draw the most common pairs on top
    sc = ax.scatter(pairs[order, 0], pairs[order, 1], c=counts[order], s=150, marker='s', cmap=cmap,
                    norm=norm, edgecolors='none', zorder=2)
    ax.set_title(title, fontsize=xs)
    ax.set_ylabel(ylabel, fontsize=xs - 2)
    ax.text(0.97, 0.03, f"Exact on {np.mean(bound == ai):.1%}", transform=ax.transAxes,
            ha='right', va='bottom', fontsize=xs - 3)
    return sc


if __name__ == "__main__":
    smiles, ai = load_table()
    ai_upper = upper_bounds(smiles)
    lower = lower_bounds(smiles)
    ai_lower, proven = lower[:, 0], lower[:, 1].astype(bool)

    if np.any(ai_upper < ai):
        raise ValueError(f"Graph RePair fell below the exact assembly index on {np.sum(ai_upper < ai)} molecules")
    if np.any(ai_lower > ai):
        raise ValueError(f"The vac bound exceeded the exact assembly index on {np.sum(ai_lower > ai)} molecules")
    print(f"Graph RePair matches the assembly index on {np.mean(ai_upper == ai):.1%} of {len(ai)} molecules, "
          f"{np.mean(ai_upper - ai):.3f} joins above it on average")
    print(f"vac matches the assembly index on {np.mean(ai_lower == ai):.1%} of {len(ai)} molecules, "
          f"{np.mean(ai - ai_lower):.3f} joins below it on average, "
          f"and proved the chain minimal for {proven.mean():.1%} within {TIMEOUT:g}s")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5), layout='constrained')
    xs = 16

    # Shared colour scale and axis limits so the panels compare directly
    cmap = LinearSegmentedColormap.from_list("blues", plt.cm.Blues(np.linspace(0.3, 1.0, 256)))
    max_count = max(np.unique(np.column_stack([ai, b]), axis=0, return_counts=True)[1].max()
                    for b in (ai_upper, ai_lower))
    norm = LogNorm(vmin=1, vmax=max_count)
    lo = min(ai.min(), ai_lower.min()) - 1
    hi = max(ai.max(), ai_upper.max()) + 1

    sc = bound_panel(ax1, ai, ai_upper, norm, cmap, "Graph RePair", 'Upper Bound', xs)
    bound_panel(ax2, ai, ai_lower, norm, cmap, "Vector Addition Chain", 'Lower Bound', xs)
    cbar = fig.colorbar(sc, ax=[ax1, ax2])
    cbar.set_label('Number of Molecules', fontsize=xs - 2)
    cbar.ax.tick_params(labelsize=xs - 4)
    fig.suptitle(f"Bounds on the Assembly Index of {len(ai):,} GDB-17 Molecules", fontsize=xs)

    for ax, label in zip([ax1, ax2], ["a)", "b)"]):
        ax.set_title(label, loc='left', fontsize=xs + 6, fontweight='bold')
        ax.plot([lo, hi], [lo, hi], '--', color='black', lw=1.5, label='Bound = Assembly Index', zorder=1)
        ax.set_xlabel('Assembly Index', fontsize=xs - 2)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect('equal')
        ax.legend(loc="upper left", fontsize=xs - 3)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.tick_params(axis='both', which='major', labelsize=xs - 2, direction='in', length=6, width=2)
        ax.tick_params(axis='both', which='both', top=True, right=True)

    plt.savefig('gdb_bounds.png', dpi=600)
    plt.savefig('gdb_bounds.pdf')
    plt.show()
