#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include "mesh.hpp"
#include <vector>

// Structure to return standard benchmark metrics across all solvers
struct SolverResult {
    std::vector<double> u;     // Computed numerical solution vector
    int iterations;            // Total iterations taken to converge
    double residual_error;     // Final residual error metric
    double computation_time;   // Wall-clock execution time in seconds
};

// 1. Successive Over-Relaxation (SOR / Gauss-Seidel) Solver
SolverResult solvePoissonSOR(const Mesh& mesh, double omega, double tol, int maxIter);

// 2. Steepest Descent Solver
SolverResult solveSteepestDescent(const Mesh& mesh, double tol, int maxIter);

// 3. Conjugate Gradient Solver
SolverResult solveConjugateGradient(const Mesh& mesh, double tol, int maxIter);

// 4. Geometric Multigrid (V-Cycle) Solver
SolverResult solveMultigrid(const Mesh& mesh, double tol, int maxIter);

#endif // SOLVERS_HPP