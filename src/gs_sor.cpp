#include "mesh.hpp"
#include "solvers.hpp"
#include "linalg.hpp"
#include <vector>
#include <cmath>
#include <chrono>

SolverResult solvePoissonSOR(const Mesh& mesh, double omega, double tol, int maxIter) {
    SolverResult result;

    int Nx = mesh.Nx;
    std::vector<double> u(Nx * Nx, 0.0);

    auto start = std::chrono::high_resolution_clock::now();
    int iter = 0;
    double maxDiff = 0.0;

    for (iter = 0; iter < maxIter; ++iter) {
        maxDiff = 0.0;
        for (int j = 0; j < Nx; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double sum = mesh.f[index(i, j, Nx)] + applyLaplacianPoint(u, i, j, Nx);

                double u_old = u[index(i, j, Nx)];
                double u_new = (1.0 - omega) * u_old + (omega / 4.0) * sum;
                u[index(i, j, Nx)] = u_new;

                maxDiff = std::max(maxDiff, std::abs(u_new - u_old));
            }
        }
        if (maxDiff < tol) break;
    }

    auto end = std::chrono::high_resolution_clock::now();

    result.u = u;
    result.iterations = iter + 1;
    result.residual_error = maxDiff;
    result.computation_time = std::chrono::duration<double>(end - start).count();

    return result;
}