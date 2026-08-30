#include "mesh.hpp"
#include "solvers.hpp"
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

namespace {

const double tol = 1e-6;
const int maxIter = 5000;
const double sorOmega = 1.9;
const char* solutionCsvPath = "output/solution.csv";
const char* solutionMetaCsvPath = "output/solution_meta.csv";
const char* benchmarkCsvPath = "output/benchmark.csv";

// The analytical solution u(x,y) = sin^2(pi x) sin^2(pi y) lies in [0, 1],
//so any converged numerical solution should too (within a small slack).
const double rangeSlack = 1e-3;

bool inValidRange(const std::vector<double>& u) {
    for (double v : u) {
        if (v < -rangeSlack || v > 1.0 + rangeSlack) return false;
    }
    return true;
}

/*
void writeRow(std::ofstream& out, const std::string& solver, int N, const SolverResult& r) {
    bool converged = r.residual_error <= tol;
    bool rangeOk = inValidRange(r.u);

    out << solver << ","
        << N << ","
        << r.iterations << ","
        << std::scientific << std::setprecision(6) << r.residual_error << ","
        << std::fixed << std::setprecision(3) << r.computation_time * 1000.0 << ","
        << (converged ? "yes" : "no") << ","
        << (rangeOk ? "yes" : "no")
        << "\n";
}
*/

void writeSolutionCSV(const Mesh& mesh, const SolverResult& r, const std::string& solver,
                       const std::string& path, const std::string& metaPath) {
    std::ofstream out(path);
    out << "x,y,u\n";
    for (int j = 0; j < mesh.Nx; ++j) {
        for (int i = 0; i < mesh.Nx; ++i) {
            double x = (i + 1) * mesh.h;
            double y = (j + 1) * mesh.h;
            out << x << "," << y << "," << r.u[index(i, j, mesh.Nx)] << "\n";
        }
    }

    std::ofstream meta(metaPath);
    meta << "solver,N,cells\n";
    meta << solver << "," << mesh.N << "," << (mesh.Nx * mesh.Nx) << "\n";
}

} // namespace

int main() {
    std::vector<int> N_values = {4, 8, 16, 32, 64, 128, 256, 512};

    std::ofstream benchmark(benchmarkCsvPath);
    benchmark << "solver,N,iterations,residual,time_ms,converged,in_range\n";

    for (int N : N_values) {
        Mesh mesh = createMesh(N);

        SolverResult sor = solvePoissonSOR(mesh, sorOmega, tol, maxIter);
        //writeRow(benchmark, "SOR", N, sor);

        SolverResult sd = solveSteepestDescent(mesh, tol, maxIter);
        //writeRow(benchmark, "SteepestD", N, sd);

        SolverResult cg = solveConjugateGradient(mesh, tol, maxIter);
        //writeRow(benchmark, "CG", N, cg);

        SolverResult mg = solveMultigrid(mesh, tol, maxIter);
        //writeRow(benchmark, "Multigrid", N, mg);

        if (N == N_values.back()) {
            writeSolutionCSV(mesh, sor, "SOR", solutionCsvPath, solutionMetaCsvPath);
        }
    }

    return 0;
}
