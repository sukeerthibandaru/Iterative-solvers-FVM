#include "mesh.hpp"
#include "solvers.hpp"
#include "linalg.hpp"
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>

namespace {

const int preSmooth    = 2;   // Pre-smoothing red-black Gauss-Seidel sweeps
const int postSmooth   = 2;   // Post-smoothing sweeps
const int coarseSweeps = 50;  // Sweeps used as the coarsest-grid "direct" solve

// One level of the grid hierarchy. Interior points only; the Dirichlet
// boundary values are implicit zeros.
struct Level {
    int Nx;                     // Interior points per dimension on this level
    double h;                   // Grid spacing on this level
    std::vector<double> u;      // Solution (or coarse-grid correction)
    std::vector<double> f;      // Unscaled RHS of this level's equation -Lap(u) = f
    std::vector<double> r;      // Residual r = f + Lap(u)

    Level(int Nx_, double h_)
        : Nx(Nx_), h(h_),
          u(Nx_ * Nx_, 0.0), f(Nx_ * Nx_, 0.0), r(Nx_ * Nx_, 0.0) {}
};

// Bounds-checked read: anything outside the interior is a zero Dirichlet boundary
inline double at(const std::vector<double>& v, int i, int j, int Nx) {
    if (i < 0 || i >= Nx || j < 0 || j >= Nx) return 0.0;
    return v[index(i, j, Nx)];
}

// Red-black Gauss-Seidel sweeps on one level
void smooth(Level& L, int sweeps) {
    const int Nx = L.Nx;
    const double h2 = L.h * L.h;

    for (int s = 0; s < sweeps; ++s) {
        for (int color = 0; color < 2; ++color) {
            for (int j = 0; j < Nx; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    if ((i + j) % 2 != color) continue;

                    double sum_nb = at(L.u, i - 1, j, Nx) + at(L.u, i + 1, j, Nx)
                                  + at(L.u, i, j - 1, Nx) + at(L.u, i, j + 1, Nx);

                    L.u[index(i, j, Nx)] = 0.25 * (sum_nb + h2 * L.f[index(i, j, Nx)]);
                }
            }
        }
    }
}

// Residual r = f + Lap(u), stored back into the level
void computeResidual(Level& L) {
    const int Nx = L.Nx;
    const double h2 = L.h * L.h;

    for (int j = 0; j < Nx; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double lap = (at(L.u, i - 1, j, Nx) + at(L.u, i + 1, j, Nx)
                        + at(L.u, i, j - 1, Nx) + at(L.u, i, j + 1, Nx)
                        - 4.0 * L.u[index(i, j, Nx)]) / h2;

            L.r[index(i, j, Nx)] = L.f[index(i, j, Nx)] + lap;
        }
    }
}

// Discrete L2 norm of a level-sized vector
double normL2(const std::vector<double>& v, double h) {
    return norm2(v) * h;
}

// Full-weighting restriction of the fine residual onto the coarse RHS.
// Coarse interior point I sits on fine interior point 2I + 1.
void restrictFullWeighting(const Level& fine, Level& coarse) {
    const int Nf = fine.Nx;
    const int Nc = coarse.Nx;

    for (int J = 0; J < Nc; ++J) {
        for (int I = 0; I < Nc; ++I) {
            int i = 2 * I + 1;
            int j = 2 * J + 1;

            double sum = 4.0 * at(fine.r, i, j, Nf)
                       + 2.0 * (at(fine.r, i - 1, j, Nf) + at(fine.r, i + 1, j, Nf)
                              + at(fine.r, i, j - 1, Nf) + at(fine.r, i, j + 1, Nf))
                       + 1.0 * (at(fine.r, i - 1, j - 1, Nf) + at(fine.r, i + 1, j - 1, Nf)
                              + at(fine.r, i - 1, j + 1, Nf) + at(fine.r, i + 1, j + 1, Nf));

            coarse.f[index(I, J, Nc)] = sum / 16.0;
        }
    }
}

// Bilinear prolongation of the coarse correction, added onto the fine solution
void prolongateAndAdd(const Level& coarse, Level& fine) {
    const int Nf = fine.Nx;
    const int Nc = coarse.Nx;

    for (int j = 0; j < Nf; ++j) {
        for (int i = 0; i < Nf; ++i) {
            double c;

            if (i % 2 == 1 && j % 2 == 1) {
                // Coincides with a coarse point: direct injection
                c = at(coarse.u, (i - 1) / 2, (j - 1) / 2, Nc);
            } else if (i % 2 == 0 && j % 2 == 1) {
                // Midpoint in x between two coarse points
                int J = (j - 1) / 2;
                c = 0.5 * (at(coarse.u, i / 2 - 1, J, Nc) + at(coarse.u, i / 2, J, Nc));
            } else if (i % 2 == 1 && j % 2 == 0) {
                // Midpoint in y between two coarse points
                int I = (i - 1) / 2;
                c = 0.5 * (at(coarse.u, I, j / 2 - 1, Nc) + at(coarse.u, I, j / 2, Nc));
            } else {
                // Cell centre: average of the four surrounding coarse points
                c = 0.25 * (at(coarse.u, i / 2 - 1, j / 2 - 1, Nc)
                          + at(coarse.u, i / 2,     j / 2 - 1, Nc)
                          + at(coarse.u, i / 2 - 1, j / 2,     Nc)
                          + at(coarse.u, i / 2,     j / 2,     Nc));
            }

            fine.u[index(i, j, Nf)] += c;
        }
    }
}

// Recursive V-cycle. levels[0] is the finest grid.
void vCycle(std::vector<Level>& levels, int lvl) {
    Level& L = levels[lvl];

    // Coarsest level: smooth it to (near) convergence
    if (lvl + 1 == static_cast<int>(levels.size())) {
        smooth(L, coarseSweeps);
        return;
    }

    smooth(L, preSmooth);
    computeResidual(L);

    Level& C = levels[lvl + 1];
    restrictFullWeighting(L, C);
    std::fill(C.u.begin(), C.u.end(), 0.0);   // Correction starts from zero

    vCycle(levels, lvl + 1);

    prolongateAndAdd(C, levels[lvl]);
    smooth(levels[lvl], postSmooth);
}

// Build the grid hierarchy from the fine mesh, halving Nx each level.
// Coarse point I sits on fine point 2I+1, so Nc = Nx/2 (integer division)
// is exactly right whether Nx is even or odd.
std::vector<Level> buildHierarchy(const Mesh& mesh) {
    std::vector<Level> levels;
    levels.emplace_back(mesh.Nx, mesh.h);

    // Mesh stores h^2 * f(x,y); the levels carry the unscaled RHS
    double h2 = mesh.h * mesh.h;
    for (size_t k = 0; k < mesh.f.size(); ++k) levels[0].f[k] = mesh.f[k] / h2;

    while (levels.back().Nx > 2) {
        int Nc = levels.back().Nx / 2;
        double hc = levels.back().h * 2.0;
        levels.emplace_back(Nc, hc);
    }

    return levels;
}

} // namespace

// Geometric multigrid (V-cycle) solver for -Lap(u) = f with zero Dirichlet data
SolverResult solveMultigrid(const Mesh& mesh, double tol, int maxIter) {
    SolverResult result;

    std::vector<Level> levels = buildHierarchy(mesh);
    Level& fine = levels[0];

    auto start = std::chrono::high_resolution_clock::now();

    // Initial residual (u == 0, so r == f)
    computeResidual(fine);
    double r0 = normL2(fine.r, fine.h);
    double res = r0;
    int iterations = 0;

    for (int cycle = 0; cycle < maxIter; ++cycle) {
        vCycle(levels, 0);
        ++iterations;

        computeResidual(fine);
        res = normL2(fine.r, fine.h);

        if (res <= tol || (r0 > 0.0 && res <= tol * r0)) break;
    }

    auto end = std::chrono::high_resolution_clock::now();

    result.u = fine.u;
    result.iterations = iterations;
    result.residual_error = res;
    result.computation_time = std::chrono::duration<double>(end - start).count();

    return result;
}
