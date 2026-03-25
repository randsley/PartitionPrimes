# Papers

Three-paper package on partition-theoretic prime detection,
extending the framework of Craig, van Ittersum & Ono (Proc. Natl. Acad. Sci. 2024).

## Papers

### Paper 1: Computational Study (main paper)
Comprehensive computational extension of the Craig-van Ittersum-Ono framework. Reproduces E1-E4, derives E5 (at d=2 in {M1,...,M5}), establishes the weight-12 barrier, recovers E6 via tau-cancellation, and confirms bounded completeness through a_max=12.

### Paper 2: The Weight-12 Cusp Obstruction
Proves that the Ramanujan tau function and the Bernoulli prime 691 force a structural obstruction at weight 12: no non-trivial prime-vanishing expression in the basis {M1, ..., M6} can involve M6. The obstruction space has dimension d+1 at each polynomial degree d.

### Paper 3: Cusp Cancellation and Recovery of E6
Shows that adjoining M7 makes the weight-12 obstruction cancellable, producing E6 as the first new prime-vanishing expression beyond E1-E5. E6 appears at polynomial degree d=4 and requires the tau-cancellation mechanism between M6 and M7.

## Companion Notebooks

Each paper has a computational companion notebook in `../notebooks/`:

| Paper | Notebook |
|-------|----------|
| Paper 1: Computational Study | `paper3_computational.ipynb` |
| Paper 2: Obstruction | `paper1_obstruction.ipynb` |
| Paper 3: Recovery | `paper2_recovery.ipynb` |

Note: The notebook numbering reflects their original creation order and does not match
the paper numbering. Each notebook's content corresponds to its named paper regardless.
