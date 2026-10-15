"""
Write the AssemblyCPP timeout lookup table read by peptide_tradeoff.py.

Runs AssemblyCPP, through assemblytheorytools, on the example's random peptides at
each timeout, and records the upper bound found and the wall-clock time of each
run. Once a run finishes exactly, or runs out of memory, longer timeouts are
skipped since they return the same result. The table is appended to as it is
written, so an interrupted run resumes where it stopped.

AssemblyCPP's memory grows quickly with peptide length, so its address space is
capped at MEMORY_LIMIT. From about 12 residues on, a run exhausts the cap within
a minute or two.

Usage: python make_peptide_tradeoff_table.py
"""
import contextlib
import csv
import io
import os
import re
import resource
import time

import assemblytheorytools as att

from peptide_tradeoff import LENGTHS, N_SAMPLES, TABLE_FILE, TIMEOUTS, heavy_atom_graph, random_peptides

FIELDS = ["length", "sample", "smiles", "timeout_s", "assembly_index", "time_s", "exact", "out_of_memory"]
# Address-space limit in bytes for this process and the AssemblyCPP runs it starts
MEMORY_LIMIT = 8_000_000_000


def assembly_index(graph, timeout: float) -> tuple:
    """
    Calculate the assembly index of a molecular graph with AssemblyCPP.

    Parameters:
        graph (nx.Graph): The heavy-atom molecular graph.
        timeout (float): Wall-clock limit in seconds.

    Returns:
        tuple: The assembly index (an upper bound if inexact, -1 if none was found),
               the wall-clock time in seconds, whether the search finished, and
               whether it ran out of memory.
    """
    log = io.StringIO()
    options = att.AssemblyCppOptions(pathway=False)
    start = time.perf_counter()
    try:
        with contextlib.redirect_stdout(log):
            ai, _, _ = att.calculate_assembly_index(graph, timeout=timeout, cpp_options=options)
    except OSError as error:
        if "bad_alloc" not in str(error):
            raise
        # The error carries the end of AssemblyCPP's log, with the best bounds found
        bounds = [int(b) for b in re.findall(r"Best assembly index: (\d+)", str(error))]
        return min(bounds, default=-1), time.perf_counter() - start, False, True
    elapsed = time.perf_counter() - start
    # assemblytheorytools reports a timeout or early stop on stdout
    exact = "timed out" not in log.getvalue() and "stopped early" not in log.getvalue()
    return ai, elapsed, exact, False


if __name__ == "__main__":
    resource.setrlimit(resource.RLIMIT_AS, (MEMORY_LIMIT, MEMORY_LIMIT))

    # Each (length, sample, timeout) already in the table, and whether it ended the sweep
    done = {}
    if os.path.isfile(TABLE_FILE):
        with open(TABLE_FILE, newline="") as f:
            done = {(int(row["length"]), int(row["sample"]), float(row["timeout_s"])):
                    row["exact"] == "True" or row["out_of_memory"] == "True"
                    for row in csv.DictReader(f)}

    with open(TABLE_FILE, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if not done:
            writer.writeheader()
        for length in LENGTHS:
            for sample, smiles in enumerate(random_peptides(length, N_SAMPLES)):
                graph = heavy_atom_graph(smiles)
                for timeout in TIMEOUTS:
                    if (length, sample, timeout) in done:
                        final = done[length, sample, timeout]
                    else:
                        ai, elapsed, exact, oom = assembly_index(graph, timeout)
                        writer.writerow(dict(length=length, sample=sample, smiles=smiles, timeout_s=f"{timeout:g}",
                                             assembly_index=ai, time_s=f"{elapsed:.6g}", exact=exact,
                                             out_of_memory=oom))
                        f.flush()
                        print(f"length {length} sample {sample} timeout {timeout:g} s: "
                              f"AI {'=' if exact else '<='} {ai} in {elapsed:.3g} s"
                              f"{' (out of memory)' if oom else ''}", flush=True)
                        final = exact or oom
                    if final:
                        break
