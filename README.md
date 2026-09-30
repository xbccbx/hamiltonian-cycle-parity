# Directed Hamiltonian-Cycle Parity

Research code and reproducibility artifacts for **Directed Hamiltonian-Cycle
Parity in O*((3/2)^n) Deterministic Time and Polynomial Space**, by Hanqing Li,
Boxuan Chen, and Xupeng Miao (Peking University).

This repository contains C++ implementations, an independent Python checker,
and archived experimental records from the authors' `v4` artifact package.
The programs compute the **parity of the number of directed Hamiltonian cycles**:
`answer = 0` means an even count and `answer = 1` means an odd count. A zero
answer does not imply that no Hamiltonian cycle exists.

## Contents

| Path | Contents |
| --- | --- |
| `artifacts/theory/code/` | Single-threaded deterministic polynomial-space implementation (`kw-ce`), Björklund–Husfeldt baseline (`bh-ce`), packed subset-DP reference, and graph generators |
| `artifacts/practical32/code/` | Parallel Las Vegas implementations, exact deduplication, and the specialized BH bipartite algorithm |
| `artifacts/verify_math.py` | Independent finite mathematical checks using only the Python standard library |
| `artifacts/build_artifacts.py` | Input-hash checks and regeneration of tables and plots from archived measurements |
| `artifacts/theory/` and `artifacts/practical32/` | Raw data, per-instance results, summaries, figures, and experimental provenance |
| `artifacts/README.md` | Detailed artifact documentation in Chinese |

## Requirements

- **C++ executables:** Linux, a C++17 compiler with OpenMP and 128-bit integer
  support (GCC is recommended), and GNU Make.
- **Mathematical checker:** Python **3.10 or newer**; no third-party packages.
- **Regenerating plots:** Python 3.10 or newer and Matplotlib; see
  `requirements.txt`.

The C++ benchmark interface uses Linux process-memory conventions. Its
`rss_kib` field is not portable to macOS without adjusting the units.

## Quick start

From the repository root:

```bash
# Independent mathematical verification (no compiler or packages needed).
python3 -B artifacts/verify_math.py

# Build the three C++ executables.
make -j2

# Run the deterministic algorithm and two references on the same small graph.
./build/bench_poly kw-ce 12 9001 dense half
./build/bench_poly bh-ce 12 9001 dense half
./build/bench_dp dp-bit 12 9001 dense 1 two

# Compare the parallel implementations on another small graph.
./build/extended kw-practical 12 1001 dense 4 10001
./build/extended bh-practical 12 1001 dense 4 10001

# Run the C++ correctness suites.
make check
```

Each benchmark invocation emits a JSON record containing `answer`, graph and
algorithm parameters, timing, and memory measurements. The examples use small
inputs so they can be tried without large memory allocations.

## Command-line interfaces

```text
bench_poly METHOD N GRAPH_SEED KIND LAYOUT
bench_dp   METHOD N GRAPH_SEED KIND THREADS [LAYOUT]
extended   METHOD N GRAPH_SEED KIND THREADS ALGORITHM_SEED
```

| Executable | Main methods | Supported vertex range |
| --- | --- | --- |
| `bench_poly` | `kw-ce`, `bh-ce` | 2–52 |
| `bench_dp` | `dp-bit` (plus research variants in `baseline.cpp`) | 2–40 |
| `extended` | `kw-practical`, `bh-practical` | 2–52 |

For `bench_poly`, `LAYOUT` is `half` or `wide`; use `half` for the archived main
comparison. It always uses one thread. For `extended`, `THREADS` is 1–32.
The supported numeric ranges are input limits, not promises that every input
fits available time or memory.

The main graph generators are `dense`, `tournament`, `bipartite`, `out2`,
`pNNN` (edge probability `NNN/1000`), and `plantedNNN` (the same probability
with a planted directed cycle). `bench_dp` uses the baseline generator; use
`dense` for the shared three-way example above. Graphs are generated from
parameters; these executables do not expose an edge-list input option.

The practical solver's algorithm seed is separate from the graph seed.
Randomness affects the computation, but the returned parity is exact.

## Reproduce the archived tables and figures

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -B artifacts/build_artifacts.py
```

The script verifies all **34** input hashes in
`artifacts/inputs_manifest.json`, checks paired records and reference
summaries, and regenerates CSV/JSON tables, LaTeX tables, and PDF/SVG/PNG plots.
It writes `artifacts/validation.json`. **It does not run performance
benchmarks.**

The mathematical checker writes `artifacts/mathematical_validation.json` and
covers 4,164 exhaustive loopless digraphs, 18 additional graphs, 5,421 ownership
graph/state pairs, and 720 translations across 24 affine-product systems.
These finite checks support correctness testing; they do not replace the
paper's proofs.

## Interpreting the experiments

The archived timings come from an Intel Core i9-14900K Linux machine. They are
historical measurements, not results measured on the reader's computer or in
continuous integration.

- **Theory configuration:** both main algorithms are single-threaded and
  deterministic with polynomial working space. The timed implementation uses
  bounded local caches, with a conservative `O(n³ log n)`-bit space bound;
  the paper's strictly streamed `O(n²)`-bit implementation is a distinct
  space claim. Packed DP is an exponential-space reference.
- **Practical configuration:** both methods use 32 hardware threads and Las
  Vegas randomization. The KW implementation uses exponential space for exact
  deduplication; the BH implementation retains polynomial space. The archived
  50-vertex KW runs reach roughly 38.2 GiB peak memory.
- The practical comparison includes 138 paired instances over 120 distinct
  graphs. Most paired timings were collected in different measurement
  periods; they are not a fully interleaved rerun. CPU frequency was not
  locked. Nine preprocessing shortcuts are excluded from the core speed
  comparison.
- Larger instances compare implementations sharing the P3 weighting routine;
  agreement there is not an independent large-instance DP verification.

See `artifacts/theory/RESULTS.md`, `artifacts/practical32/RESULTS.md`, and their
`IMPLEMENTATION_NOTES.md` files for the full qualifications and per-instance
results. Historical source paths and build hashes in the provenance records
describe the original experiments; the archived remote binaries and original
benchmark orchestration scripts are not included.

## Source provenance and citation

The C++ sources, mathematical checker, plot-generation script, and frozen
inputs are preserved from `v4.zip`. Repository documentation and the build
and CI entry points were added for this public release. The original package's
manuscript-delivery report is omitted because it refers to manuscript files
outside this code release.

The BH reference is implemented in the included benchmark source; attribution
to the Björklund–Husfeldt algorithms is retained in the code and implementation
notes. This repository does not bundle an external BH software distribution.

Please cite the accompanying manuscript when using this work. Author and
title metadata are provided in `CITATION.cff`; no publication venue or DOI is
asserted in this repository.

## License

MIT License. See `LICENSE`.

## 中文说明

本仓库公开论文的 C++ 实现、Python 数学验证、实验数据及图表复算脚本。
论文作者为 Hanqing Li、Boxuan Chen、Xupeng Miao。详细中文说明见
[`artifacts/README.md`](artifacts/README.md)。

先用 Python 3.10 或更新版本运行 `python3 -B artifacts/verify_math.py`；
在带有 GCC/OpenMP 的 Linux 环境中运行 `make` 编译、`make check` 检查。
`answer` 表示 Hamiltonian 环数量的奇偶性，输出 0 不代表不存在 Hamiltonian 环。
已有性能数字来自历史实验；复算表图不会重新测量性能。
