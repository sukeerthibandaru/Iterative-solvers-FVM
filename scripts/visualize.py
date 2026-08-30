"""Plot the numerical Poisson solution against the analytical solution.

Reads output/solution.csv (written by main.cpp: columns x,y,u) and
output/solution_meta.csv (columns solver,N,cells), and plots the numerical
solution, the analytical solution u(x,y) = sin^2(pi x) sin^2(pi y), and
their difference as heatmaps.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DATA_PATH = Path(__file__).resolve().parent.parent / "output" / "solution.csv"
META_PATH = Path(__file__).resolve().parent.parent / "output" / "solution_meta.csv"
OUT_PATH = Path(__file__).resolve().parent.parent / "output" / "heatmap_comparison.png"


def load_grid(path):
    df = pd.read_csv(path)
    pivot = df.pivot(index="y", columns="x", values="u")
    x = pivot.columns.to_numpy()
    y = pivot.index.to_numpy()
    u_num = pivot.to_numpy()
    return x, y, u_num


def load_meta(path):
    row = pd.read_csv(path).iloc[0]
    return str(row["solver"]), int(row["N"]), int(row["cells"])


def analytical(x, y):
    X, Y = np.meshgrid(x, y)
    return np.sin(np.pi * X) ** 2 * np.sin(np.pi * Y) ** 2


def main():
    x, y, u_num = load_grid(DATA_PATH)
    solver, N, cells = load_meta(META_PATH)
    u_exact = analytical(x, y)
    error = u_num - u_exact

    l2_err = np.sqrt(np.mean(error ** 2))
    max_err = np.max(np.abs(error))
    print(f"Grid: {u_num.shape[1]} x {u_num.shape[0]} interior points")
    print(f"L2 error  : {l2_err:.6e}")
    print(f"Max error : {max_err:.6e}")

    extent = [x.min(), x.max(), y.min(), y.max()]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    im0 = axes[0].imshow(u_num, origin="lower", extent=extent, cmap="viridis")
    axes[0].set_title("Numerical solution")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(u_exact, origin="lower", extent=extent, cmap="viridis")
    axes[1].set_title("Analytical solution")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    err_bound = max(max_err, 1e-15)
    im2 = axes[2].imshow(error, origin="lower", extent=extent, cmap="coolwarm",
                          vmin=-err_bound, vmax=err_bound)
    axes[2].set_title("Error (numerical - analytical)")
    axes[2].set_xlabel("x")
    axes[2].set_ylabel("y")
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    fig.suptitle(
        r"Poisson solver vs analytical $u(x,y)=\sin^2(\pi x)\sin^2(\pi y)$"
        f"\nSolver: {solver}   N: {N}   Cells: {cells}"
    )
    fig.tight_layout()

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, dpi=150)
    print(f"Saved heatmap comparison to {OUT_PATH}")


if __name__ == "__main__":
    main()
