#include "mesh.hpp"
#include "solvers.hpp"
#include "linalg.hpp"
#include <vector>
#include <cmath>
#include <chrono>

// Conjugate gradient solver for -Lap(u) = f with zero Dirichlet data
SolverResult solveConjugateGradient(const Mesh& mesh, double tol, int maxIter) {
    SolverResult result;

    const int Nx = mesh.Nx;
    const int n = Nx * Nx;

    std::vector<double> x(n, 0.0);
    std::vector<double> r = mesh.f;   // r = b - A*x, x0 = 0 so r0 = b
    std::vector<double> p = r;
    std::vector<double> Ap(n);

    auto start = std::chrono::high_resolution_clock::now();

    double resNorm = norm2(r);
    double resNorm0 = resNorm;
    int iter = 0;

    while (resNorm > tol * resNorm0 && iter < maxIter) {
        applyLaplacian(p, Ap, Nx);

        double rr = dot(r, r);
        double pAp = dot(p, Ap);
        if (std::fabs(pAp) < 1e-12) break;
        double alpha = rr / pAp;

        for (int k = 0; k < n; ++k) {
            x[k] += alpha * p[k];
            r[k] -= alpha * Ap[k];
        }

        double rrNew = dot(r, r);
        double beta = rrNew / rr;

        for (int k = 0; k < n; ++k) {
            p[k] = r[k] + beta * p[k];
        }

        resNorm = std::sqrt(rrNew);
        ++iter;
    }

    auto end = std::chrono::high_resolution_clock::now();

    result.u = x;
    result.iterations = iter;
    result.residual_error = resNorm;
    result.computation_time = std::chrono::duration<double>(end - start).count();

    return result;
}
