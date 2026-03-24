# Papers

LaTeX sources for the three-paper package on partition-theoretic prime detection,
extending the framework of Craig, van Ittersum & Ono (Proc. Natl. Acad. Sci. 2024).

## Papers

### Paper 1: The Weight-12 Cusp Obstruction
**File:** `paper1-obstruction.tex`

Proves that the Ramanujan tau function and the Bernoulli prime 691 force a structural obstruction at weight 12: no non-trivial prime-vanishing expression in the basis {M1, ..., M6} can involve M6. The obstruction space has dimension d+1 at each polynomial degree d.

### Paper 2: Cusp Cancellation and Recovery of E6
**File:** `paper2-recovery.tex`

Shows that adjoining M7 makes the weight-12 obstruction cancellable, producing E6 as the first new prime-vanishing expression beyond E1-E5. E6 appears at polynomial degree d=4 and requires the tau-cancellation mechanism between M6 and M7.

### Paper 3: Computational Study (main paper)
**File:** `paper.tex`

Comprehensive computational extension of the Craig-van Ittersum-Ono framework. Reproduces E1-E4, derives E5 (at d=2 in {M1,...,M5}), establishes the weight-12 barrier, recovers E6 via tau-cancellation, and confirms bounded completeness through a_max=12.

### Supporting Files

| File | Description |
|------|-------------|
| `appendix-computations.tex` | Appendix A: explicit formulas for E5 and E6, extended search tables |
| `shared-preamble.tex` | Common LaTeX preamble shared across all three papers |
| `three_paper_insertions.tex` | Cross-reference insertions linking the papers |

## Building

The papers use standard LaTeX with no special build system. Compile with:

```bash
cd Paper
pdflatex paper.tex              # Main paper (may need 2 runs for references)
pdflatex paper1-obstruction.tex
pdflatex paper2-recovery.tex
```

All three papers share `shared-preamble.tex` via `\input{}`.

## Companion Notebooks

Each paper has a computational companion notebook in `../notebooks/`:

| Paper | Notebook |
|-------|----------|
| Paper 1 | `paper1_obstruction.ipynb` |
| Paper 2 | `paper2_recovery.ipynb` |
| Paper 3 | `paper3_computational.ipynb` |
