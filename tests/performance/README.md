# Performance tests

Benchmarks for matid's structure-based clustering (`SBC`) that measure CPU
time and peak memory usage across a range of system sizes for both ordered
(crystalline FCC copper) and disordered (random copper) inputs. Results are
tracked per commit so individual commits can be benchmarked and compared.

## What is measured

For each system size and category, a single classification is run, namely
`SBC().get_clusters(atoms)`.

- **CPU time** — wall-clock seconds around the call.
- **Peak memory** — `ru_maxrss` minus the resident-set size at the start of
  the call, in MiB.

The test systems:

- **Ordered** — bulk FCC copper (a = 3.6 Å, cubic 4-atom unit cell) tiled to
  approximate the requested size. The actual atom count is rounded down to
  the nearest perfect cube of unit cells (e.g. requesting 1500 yields 1372
  atoms).
- **Unordered** — copper atoms placed at uniformly random scaled positions in
  a cubic cell sized to match the ordered density. The atom count is exact.

A fixed NumPy seed makes the unordered systems reproducible across runs.

## Output layout

Results are written to a per-commit folder:

```
results/<version_id>/
    results.json    # raw timings and memory measurements
    results.pdf     # plot rendered by plot.sh
```

`<version_id>` is chosen as follows:

| Situation                          | Folder name                   |
|------------------------------------|-------------------------------|
| Inside the matid git tree, clean   | git short SHA, e.g. `24a8249` |
| Inside the matid git tree, dirty   | `<sha>-dirty`                 |
| matid installed without a git tree | `v<version>`, e.g. `v2.1.6`   |

The `-dirty` suffix prevents a benchmark run with uncommitted changes from
silently overwriting a clean-commit result.

## Running

From this directory, with matid installed in the active Python environment:

```sh
./benchmark_cpu.sh       # CPU timings for ordered and unordered systems
./benchmark_memory.sh    # Peak memory for ordered and unordered systems
./plot.sh                # Renders results.pdf from results.json
```

Each shell script loops over a fixed set of sizes (ordered: 500–4000 step
500; unordered: 100–500 step 100) and invokes `benchmark-cpu-single` or
`benchmark-memory-single` once per size. Sizes already present in
`results.json` are skipped, so runs are idempotent and a partial run can be
resumed by re-invoking the script.

`plot.sh` tolerates partial data: if only CPU or only memory measurements
are present, the corresponding panel is plotted and the other is left empty.
If no data exists for the current `version_id`, it exits with a message
telling you to run the benchmark scripts first.

The shell scripts call `python` directly. If matid lives in a virtual
environment that is not on `PATH`, activate it first:

```sh
source ../../.venv/bin/activate
```

## Benchmarking a specific commit

```sh
git checkout <commit>
uv sync                      # or: pip install -e .
cd tests/performance
./benchmark_cpu.sh
./benchmark_memory.sh
./plot.sh
# → results/<short-sha>/results.{json,pdf}
```

Repeat for each commit to build up a comparison set. Plot titles include
both the matid version and the SHA, so generated PDFs are self-identifying.

## Adding sizes

To extend the benchmarks, append values to the loops in `benchmark_cpu.sh`
and `benchmark_memory.sh`. Existing measurements in `results.json` are
preserved.
