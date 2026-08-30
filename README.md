# Iterative-solvers-FVM

Finite-volume solvers for the 2D Poisson equation on the unit square, benchmarked
against a manufactured analytical solution.

## Problem

Solves

```
-Δu(x, y) = f(x, y),   (x, y) ∈ (0, 1)²
u = 0 on the boundary
```

with `f` chosen so the exact solution is

```
u(x, y) = sin²(πx) · sin²(πy)
```

The domain is discretized with a uniform 5-point-stencil finite-volume grid
(`Mesh`, built by [createMesh](src/mesh.cpp)); `Nx = N - 2` interior points
per dimension, spacing `h = 1 / (N - 1)`.

## Solvers

All four solvers share the `SolverResult { u, iterations, residual_error, computation_time }`
interface declared in [solvers.hpp](include/solvers.hpp):

| # | Method | File |
|---|--------|------|
| 1 | Successive Over-Relaxation (SOR / Gauss-Seidel) | [gs_sor.cpp](src/gs_sor.cpp) |
| 2 | Steepest Descent | [steepest_descent.cpp](src/steepest_descent.cpp) |
| 3 | Conjugate Gradient | [conjugate_gradient.cpp](src/conjugate_gradient.cpp) |
| 4 | Geometric Multigrid (V-cycle, red-black GS smoothing) | [multigrid.cpp](src/multigrid.cpp) |

Shared discrete-Laplacian operators (matrix-free `A·v`, dot product, norm) live in
[linalg.hpp](include/linalg.hpp) / [linalg.cpp](src/linalg.cpp).

## Project layout

```
include/        Mesh and solver headers
src/            Solver implementations + main.cpp driver
scripts/        visualize.py — plots numerical vs. analytical solution
sample_results/ Example output (CSV + heatmap) from a previous run
docs/report.pdf Write-up
```

## Building

There is no build system yet — compile directly with any C++17 compiler. To build
the benchmark driver ([main.cpp](src/main.cpp)), which runs every solver over a
sweep of grid sizes and writes results to `output/`:

```bash
g++ -std=c++17 -O2 -Iinclude src/*.cpp -o poisson_solver
mkdir -p output
./poisson_solver
```

This writes `output/benchmark.csv` (timings/iterations per solver and grid size)
and, for the finest `N`, `output/solution.csv` + `output/solution_meta.csv` (the
SOR solution field, for plotting).

## Visualizing a solution

After running the driver (or any program that writes `output/solution.csv` +
`output/solution_meta.csv` in the same format), plot the numerical solution
against the analytical one:

```bash
pip install numpy pandas matplotlib
python scripts/visualize.py
```

This saves a three-panel heatmap (numerical, analytical, error) to
`output/heatmap_comparison.png` and prints the L2/max error. An example is in
[sample_results/heatmap_comparison.png](sample_results/heatmap_comparison.png).

## Notes

- SOR/Gauss-Seidel converges slowly at high resolution; the default
  `maxIter = 5000` in `main.cpp` may not be enough to reach `tol = 1e-6` for
  large `N` (e.g. 512), in which case the residual plateaus above tolerance.
  Conjugate Gradient and Multigrid converge in far fewer iterations at scale.
- `sorOmega = 1.9` in `main.cpp` is tuned for this problem; the optimal
  relaxation factor depends on grid size.
