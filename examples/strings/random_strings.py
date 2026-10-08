import csv
import gc
import os
import time

import matplotlib.pyplot as plt
import numpy as np

import assemblycfg as CFG

# Setting plot aesthetics for better visibility
plt.rcParams['axes.linewidth'] = 2.0

ALPHABET = list("abcd")
# String lengths for the assembly index panel and the compute time panel
AI_LENGTHS = np.linspace(10, 100, 20).astype(int)
TIME_LENGTHS = np.unique(np.geomspace(10, 1000, 50).astype(int))
N_SAMPLES = 20
SEED = 2024
# RePair and the LZ bound are timed as the fastest of several calls on each string
REPEATS = 5
# Exact assembly indices and compute times from AssemblyCPP, written by an offline script
TABLE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "random_strings_ai.csv")


def random_strings(length: int, n_samples: int = N_SAMPLES, seed: int = SEED) -> list:
    """
    Sample random strings uniformly from a 4-character alphabet.

    The strings depend only on the length and seed, so the lookup table and this
    example draw the same strings.

    Parameters:
        length (int): The length of each string.
        n_samples (int): The number of strings to sample.
        seed (int): The random seed.

    Returns:
        list of str: The sampled strings.
    """
    rng = np.random.default_rng([seed, length])
    return ["".join(s) for s in rng.choice(ALPHABET, size=(n_samples, length))]


def timed(func, s):
    """
    Call a function on a string and time it.

    Garbage collection is paused during the call, as in timeit, so that a collection
    of the whole interpreter's objects is not charged to whichever call triggers it.

    Parameters:
        func (callable): The function to call.
        s (str): The input string.

    Returns:
        tuple: The function's result and the wall-clock time in seconds.
    """
    gc_was_enabled = gc.isenabled()
    gc.disable()
    try:
        start = time.perf_counter()
        result = func(s)
        elapsed = time.perf_counter() - start
    finally:
        if gc_was_enabled:
            gc.enable()
    return result, elapsed


def best_of(func, s, repeats=REPEATS):
    """
    Call a function on a string several times and keep the fastest time.

    Parameters:
        func (callable): The function to call.
        s (str): The input string.
        repeats (int): The number of calls.

    Returns:
        tuple: The function's result and the fastest wall-clock time in seconds.
    """
    runs = [timed(func, s) for _ in range(repeats)]
    return runs[0][0], min(t for _, t in runs)


def repair_upper_bound(s: str) -> int:
    """Upper bound on the assembly index from RePair."""
    return CFG.repair_with_pathways(s)[0]


def load_table(file: str = TABLE_FILE) -> dict:
    """
    Load the AssemblyCPP lookup table.

    Parameters:
        file (str): Path to the CSV lookup table.

    Returns:
        tuple: A dict mapping each string length to a list of (string, assembly index, time, exact)
               rows, and the AssemblyCPP timeout in seconds. The assembly index is an upper bound
               when the calculation timed out.
    """
    table = {}
    timeouts = set()
    with open(file, newline="") as f:
        for row in csv.DictReader(f):
            entry = (row["string"], int(row["assembly_index"]), float(row["time_s"]), row["exact"] == "True")
            table.setdefault(int(row["length"]), []).append(entry)
            timeouts.add(float(row["timeout_s"]))
    if len(timeouts) != 1:
        raise ValueError(f"The lookup table mixes AssemblyCPP timeouts: {sorted(timeouts)}")
    return table, timeouts.pop()


def band(ax, x, samples, color, style, label, lw):
    """Plot the mean of each row of samples with a min-max band."""
    samples = np.asarray(samples)
    ax.fill_between(x, samples.min(axis=1), samples.max(axis=1), color=color, alpha=0.25, lw=0)
    ax.plot(x, samples.mean(axis=1), style, color=color, lw=lw, label=label)


if __name__ == "__main__":
    table, timeout = load_table()

    # Bounds on every sample, and AssemblyCPP on the samples in the table (fewer at long
    # lengths, where it always times out)
    repair_ai, lz_ai, acpp_ai, acpp_exact = [], [], [], []
    for length in AI_LENGTHS:
        strings = random_strings(length)
        rows = table[length]
        repair_ai.append(np.mean([repair_upper_bound(s) for s in strings]))
        lz_ai.append(np.mean([CFG.lz_lower_bound(s) for s in strings]))
        acpp_ai.append(np.mean([ai for _, ai, _, _ in rows]))
        acpp_exact.append(all(exact for *_, exact in rows))
    acpp_ai = np.array(acpp_ai)
    acpp_exact = np.array(acpp_exact)

    # Compute times of the bounds out to long strings, as the fastest of several sweeps over
    # every length. Spreading the repeats over the run, rather than repeating each call back
    # to back, keeps a slow stretch of the machine from raising every repeat of one length.
    strings = [random_strings(length) for length in TIME_LENGTHS]
    repair_time = np.full((len(TIME_LENGTHS), N_SAMPLES), np.inf)
    lz_time = np.full((len(TIME_LENGTHS), N_SAMPLES), np.inf)
    for _ in range(REPEATS):
        for i, samples in enumerate(strings):
            for j, s in enumerate(samples):
                repair_time[i, j] = min(repair_time[i, j], timed(repair_upper_bound, s)[1])
                lz_time[i, j] = min(lz_time[i, j], timed(CFG.lz_lower_bound, s)[1])

    # AssemblyCPP compute times wherever every sample finished exactly
    acpp_lengths = [length for length in TIME_LENGTHS
                    if length in table and all(exact for *_, exact in table[length])]
    acpp_time = [[t for _, _, t, _ in table[length]] for length in acpp_lengths]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    xs = 16

    ax1.plot(AI_LENGTHS, repair_ai, '--', color='black', lw=1.5, label='RePair Upper Bound')
    ax1.plot(AI_LENGTHS[~acpp_exact], acpp_ai[~acpp_exact], 'v', color='red', ms=8, ls='none',
             label=f'AssemblyCPP {timeout:g}s Timeout Upper Bound')
    ax1.plot(AI_LENGTHS[acpp_exact], acpp_ai[acpp_exact], 'P', color='green', ms=9, ls='none',
             label='AssemblyCPP Exact Assembly Index')
    ax1.plot(AI_LENGTHS, lz_ai, ':', color='black', lw=1.5, label='LZ Lower Bound')
    ax1.set_title("Assembly Index vs String Length", fontsize=xs)
    ax1.set_xlabel('String Length', fontsize=xs - 2)
    ax1.set_ylabel('Joining Operations', fontsize=xs - 2)
    ax1.set_xlim(AI_LENGTHS[0], AI_LENGTHS[-1])
    ax1.legend(loc="upper left", fontsize=xs - 3)

    band(ax2, TIME_LENGTHS, repair_time, 'black', '--', 'RePair Upper Time', lw=3)
    band(ax2, acpp_lengths, acpp_time, 'green', '-', 'AssemblyCPP Exact Time', lw=3)
    band(ax2, TIME_LENGTHS, lz_time, 'black', ':', 'LZ Lower Time', lw=3)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_title("Compute Time vs String Length", fontsize=xs)
    ax2.set_xlabel('String Length', fontsize=xs - 2)
    ax2.set_ylabel('Compute Time (s)', fontsize=xs - 2)
    ax2.set_xlim(TIME_LENGTHS[0], TIME_LENGTHS[-1])
    # Crop at the slowest mean AssemblyCPP time
    ax2.set_ylim(top=np.mean(acpp_time, axis=1).max())
    ax2.legend(loc="upper right", fontsize=xs - 3)

    for ax in [ax1, ax2]:
        ax.tick_params(axis='both', which='major', labelsize=xs - 2, direction='in', length=6, width=2)
        ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1)
        ax.tick_params(axis='both', which='both', top=True, right=True)

    plt.tight_layout()
    plt.savefig('random_strings.png', dpi=600)
    plt.savefig('random_strings.pdf')
    plt.show()
