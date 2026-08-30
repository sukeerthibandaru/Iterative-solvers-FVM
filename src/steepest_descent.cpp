#include "mesh.hpp"
#include "solvers.hpp"
#include "linalg.hpp"
#include <vector>
#include <cmath>
#include <chrono>

// Steepest descent solver for -Lap(u) = f with zero Dirichlet data
SolverResult solveSteepestDescent(const Mesh& mesh, double tol, int maxIter) {
    SolverResult result;

    const int Nx = mesh.Nx;
    const int n = Nx * Nx;

    std::vector<double> x(n, 0.0);
    std::vector<double> r = mesh.f;   // r = b - A*x, x0 = 0 so r0 = b
    std::vector<double> Ar(n);

    auto start = std::chrono::high_resolution_clock::now();

    double resNorm = norm2(r);
    int iter = 0;

    while (resNorm > tol && iter < maxIter) {
        applyLaplacian(r, Ar, Nx);
        double rr = dot(r, r);
        double rAr = dot(r, Ar);
        if (std::fabs(rAr) < 1e-12) break;
        double alpha = rr / rAr;

        for (int k = 0; k < n; ++k) {
            x[k] += alpha * r[k];
            r[k] -= alpha * Ar[k];
        }

        resNorm = norm2(r);
        ++iter;
    }

    auto end = std::chrono::high_resolution_clock::now();

    result.u = x;
    result.iterations = iter;
    result.residual_error = resNorm;
    result.computation_time = std::chrono::duration<double>(end - start).count();

    return result;
}
