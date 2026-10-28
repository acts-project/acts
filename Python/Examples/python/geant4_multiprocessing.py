"""Run independent Geant4 event ranges in fresh processes.

Call from a ``if __name__ == "__main__"`` guard in a script. Workers must be
importable functions (or picklable partials) that construct their own detector,
random-number service and single-threaded sequencer. Never pass Geant4 or ACTS
geometry objects from the parent process.
"""

from multiprocessing import get_context
from multiprocessing.connection import wait
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory


def runGeant4Processes(worker, *, events, processes, skip=0):
    """Call ``worker(begin, end)`` concurrently for disjoint event ranges.

    Each spawned process executes exactly one range, keeping Geant4's lifetime
    local to that process. Use the same random seed in every worker and set the
    sequencer's ``skip=begin``, ``events=end-begin`` and ``numThreads=1``.
    Returns the ranges in event order. Any failed worker fails the whole job
    and terminates remaining workers, including on interrupts.
    """
    if events < 0 or skip < 0 or processes < 1:
        raise ValueError(
            "events and skip must be nonnegative; processes must be positive"
        )
    count = min(processes, max(events, 1))
    size, remainder = divmod(events, count)
    ranges = []
    begin = skip
    for index in range(count):
        end = begin + size + (index < remainder)
        ranges.append((begin, end))
        begin = end

    context = get_context("spawn")
    children = []
    try:
        for begin, end in ranges:
            child = context.Process(
                target=worker, args=(begin, end), name=f"Geant4-events-{begin}-{end}"
            )
            child.start()
            children.append(child)
        pending = {child.sentinel: child for child in children}
        while pending:
            for sentinel in wait(pending):
                child = pending.pop(sentinel)
                child.join()
                if child.exitcode != 0:
                    raise RuntimeError(
                        f"{child.name} failed with exit code {child.exitcode}"
                    )
    finally:
        for child in children:
            if child.is_alive():
                child.terminate()
        for child in children:
            child.join()
            child.close()
    return ranges


def mergeGeant4RootFiles(inputs, output):
    """Merge worker ROOT files in event order, publishing only on success."""
    hadd = shutil.which("hadd")
    if hadd is None:
        raise RuntimeError("ROOT's hadd must be on PATH to merge Geant4 worker output")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="geant4-merge-", dir=output.parent) as temporary:
        merged = Path(temporary) / output.name
        subprocess.run([hadd, "-f", str(merged), *map(str, inputs)], check=True)
        merged.replace(output)
