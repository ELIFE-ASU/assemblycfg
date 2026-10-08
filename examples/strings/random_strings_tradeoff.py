import csv
import os
from string import ascii_lowercase

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.container import ErrorbarContainer
from matplotlib.ticker import NullLocator, PercentFormatter

import assemblycfg as CFG
from random_strings import best_of, random_strings, repair_upper_bound

# Setting plot aesthetics for better visibility
plt.rcParams['axes.linewidth'] = 2.0

# String lengths, one panel each, in rows of PANELS_PER_ROW
LENGTHS = [50, 100, 300, 1000, 3000]
PANELS_PER_ROW = 3
N_SAMPLES = 10
# AssemblyCPP timeouts in seconds. A run stopped by the timeout returns the best upper bound found.
TIMEOUTS = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100]
# AssemblyCPP results at each timeout, written by an offline script
TABLE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "random_strings_tradeoff.csv")


def load_table(file: str = TABLE_FILE) -> dict:
    """
    Load the AssemblyCPP timeout lookup table.

    Runs at timeouts beyond one that finished exactly are not in the table, and take the
    result of the exact run.

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
                                             row["exact"] == "True")
    for runs in table.values():
        for timeout in TIMEOUTS:
            if timeout not in runs and any(exact for _, _, exact in runs.values()):
                runs[timeout] = next(run for run in runs.values() if run[2])
    return table


def panel_grid(fig, n_panels, per_row=PANELS_PER_ROW):
    """
    Add panels with shared axes in rows, centering a shorter last row.

    Parameters:
        fig (Figure): The figure to add the panels to.
        n_panels (int): The number of panels.
        per_row (int): The number of panels in a full row.

    Returns:
        list of tuple: Each panel's axes and whether it starts its row.
    """
    n_rows = -(-n_panels // per_row)
    # Each panel spans two grid columns, so a short row can be offset by half a panel
    grid = fig.add_gridspec(n_rows, 2 * per_row)
    panels = []
    for i in range(n_panels):
        row, col = divmod(i, per_row)
        start = 2 * col + per_row - min(per_row, n_panels - row * per_row)
        shared = panels[0][0] if panels else None
        ax = fig.add_subplot(grid[row, start:start + 2], sharex=shared, sharey=shared)
        panels.append((ax, col == 0))
    return panels


def point(ax, times, errors, marker, label, color='black'):
    """Plot the mean time and error of one method, with min-max whiskers over the samples."""
    times, errors = np.asarray(times), np.asarray(errors)
    x, y = times.mean(), errors.mean()
    ax.errorbar(x, y, xerr=[[x - times.min()], [times.max() - x]],
                yerr=[[y - errors.min()], [errors.max() - y]],
                fmt=marker, color=color, ms=11, lw=1.5, capsize=3, label=label)


if __name__ == "__main__":
    table = load_table()

    n_rows = -(-len(LENGTHS) // PANELS_PER_ROW)
    fig = plt.figure(figsize=(4.4 * PANELS_PER_ROW, 4.3 * n_rows + 0.6))
    panels = panel_grid(fig, len(LENGTHS))
    xs = 16

    for (ax, row_start), length, letter in zip(panels, LENGTHS, ascii_lowercase):
        strings = random_strings(length, N_SAMPLES)
        runs = [table[length, sample] for sample in range(N_SAMPLES)]

        repair = [best_of(repair_upper_bound, s) for s in strings]
        lz = [best_of(CFG.lz_lower_bound, s) for s in strings]

        # The exact assembly index where AssemblyCPP finished, else the best upper bound found
        n_exact = sum(any(e for _, _, e in r.values()) for r in runs)
        reference = np.array([min([ai for ai, _, _ in r.values() if ai >= 0] + [rp])
                              for r, (rp, _) in zip(runs, repair)])

        def error(bounds):
            return (np.asarray(bounds) - reference) / reference

        if n_exact < N_SAMPLES:
            # Where the true assembly index lies: at the reference when exact, else between the
            # LZ lower bound and the best upper bound
            lower = [ref if any(e for _, _, e in r.values()) else b
                     for r, ref, (b, _) in zip(runs, reference, lz)]
            ax.axhspan(error(lower).mean(), 0, color='grey', alpha=0.15, lw=0, label='Range of True Value')
        ax.axhline(0, color='grey', lw=1)

        # AssemblyCPP's best bound within each timeout. Its search is deterministic, so a bound found
        # within a shorter timeout is also found within a longer one. assemblytheorytools sometimes
        # returns no bound from a run stopped by the timeout, and before any bound is found the
        # trivial bound, length - 1, holds. Runs that return no bound overrun the timeout while
        # stopping, so a run's time is capped at its timeout.
        acpp_time = [np.mean([min(r[timeout][1], timeout) for r in runs]) for timeout in TIMEOUTS]
        acpp_err = []
        for i in range(len(TIMEOUTS)):
            best = [min([r[t][0] for t in TIMEOUTS[:i + 1] if r[t][0] >= 0] + [length - 1]) for r in runs]
            acpp_err.append(error(best))
        acpp_err = np.array(acpp_err)
        ax.fill_between(acpp_time, acpp_err.min(axis=1), acpp_err.max(axis=1), color='green', alpha=0.25, lw=0)
        ax.plot(acpp_time, acpp_err.mean(axis=1), '-o', color='green', lw=2.5, ms=6, label='AssemblyCPP')

        point(ax, [t for _, t in repair], error([b for b, _ in repair]), 's', 'RePair Upper Bound')
        point(ax, [t for _, t in lz], error([b for b, _ in lz]), 'v', 'LZ Lower Bound')

        ax.set_xscale('log')
        # Linear within 10% of the reference, where the bounds are compared, and logarithmic beyond
        ax.set_yscale('symlog', linthresh=0.1, linscale=2)
        ax.set_yticks([-0.3, -0.05, 0, 0.05, 0.1, 1])
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_title(f"Length {length} ({n_exact}/{N_SAMPLES} Exact)", fontsize=xs - 2)
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
    plt.savefig('random_strings_tradeoff.png', dpi=600)
    plt.savefig('random_strings_tradeoff.pdf')
    plt.show()
