import gc
import tracemalloc

import matplotlib.pyplot as plt
import numpy as np

import assemblycfg as CFG
from random_strings import random_strings

# Setting plot aesthetics for better visibility
plt.rcParams['axes.linewidth'] = 2.0

LENGTHS = np.unique(np.geomspace(10, 100_000, 25).astype(int))
N_SAMPLES = 10
# Exponents are fitted to the lengths of the last decade measured
FIT_DECADES = 1


def repair_bound_only(s: str) -> int:
    """Upper bound on the assembly index from RePair, without building the pathway."""
    return CFG.ai_core(s)[0]


# Label, function and line style
METHODS = [
    ('RePair Upper Bound', repair_bound_only, '--'),
    ('LZ Lower Bound', CFG.lz_lower_bound, ':'),
]


def peak_memory(func, s):
    """
    Measure the peak memory a function allocates while it runs on a string.

    Uses tracemalloc, which counts every allocation made through Python's allocators,
    so it includes the result but not memory allocated before the call, such as the
    input. The call runs several times slower while traced.

    Parameters:
        func (callable): The function to call.
        s (str): The input string.

    Returns:
        int: The peak traced memory during the call, in bytes.
    """
    gc.collect()
    tracemalloc.start()
    try:
        func(s)
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def exponent(lengths, peaks):
    """Fit peak memory to a power of the length over the last FIT_DECADES decades."""
    lengths, peaks = np.asarray(lengths), np.asarray(peaks).mean(axis=1)
    fit = lengths >= lengths[-1] / 10 ** FIT_DECADES
    return np.polyfit(np.log(lengths[fit]), np.log(peaks[fit]), 1)[0]


def band(ax, x, samples, style, label):
    """Plot the mean of each row of samples with a min-max band."""
    samples = np.asarray(samples)
    ax.fill_between(x, samples.min(axis=1), samples.max(axis=1), color='black', alpha=0.2, lw=0)
    ax.plot(x, samples.mean(axis=1), style, color='black', lw=2.5, label=label)


if __name__ == "__main__":
    strings = {length: random_strings(length, N_SAMPLES) for length in LENGTHS}

    results = []
    for label, func, style in METHODS:
        peaks = np.array([[peak_memory(func, s) for s in strings[length]] for length in LENGTHS])
        k = exponent(LENGTHS, peaks)
        print(f"{label}: {peaks[-1].mean() / 1e6:.1f} MB at length {LENGTHS[-1]}, "
              f"peak memory ~ length^{k:.2f}", flush=True)
        results.append((label, k, style, peaks))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    xs = 16

    # The exponents describe the total peak memory, so only its legend shows them
    for label, k, style, peaks in results:
        band(ax1, LENGTHS, peaks / 1e6, style, rf'{label} ($\propto n^{{{k:.2f}}}$)')
        band(ax2, LENGTHS, peaks / LENGTHS[:, None], style, label)

    ax1.set_title("Peak Memory vs String Length", fontsize=xs)
    ax1.set_ylabel('Peak Memory (MB)', fontsize=xs - 2)
    ax2.set_title("Peak Memory per Character vs String Length", fontsize=xs)
    ax2.set_ylabel('Peak Memory per Character (bytes)', fontsize=xs - 2)

    for ax in [ax1, ax2]:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(LENGTHS[0], LENGTHS[-1])
        ax.set_xlabel('String Length', fontsize=xs - 2)
        ax.legend(loc="upper left", fontsize=xs - 3)
        ax.tick_params(axis='both', which='major', labelsize=xs - 2, direction='in', length=6, width=2)
        ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1)
        ax.tick_params(axis='both', which='both', top=True, right=True)

    # Headroom above the curves for the legend, and on the per-character axis, room
    # below down to a whole decade so that at least two decades are labelled
    ax1.set_ylim(top=ax1.get_ylim()[1] * 4)
    bottom, top = ax2.get_ylim()
    ax2.set_ylim(10 ** np.floor(np.log10(bottom)), top * 2)

    plt.tight_layout()
    plt.savefig('random_strings_memory_scaling.png', dpi=600)
    plt.savefig('random_strings_memory_scaling.pdf')
    plt.show()
