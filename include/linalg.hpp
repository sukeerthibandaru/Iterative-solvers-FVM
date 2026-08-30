#ifndef LINALG_HPP
#define LINALG_HPP

#include <vector>

// Apply the discrete negative Laplacian (5-point stencil, zero Dirichlet
// boundaries) to v on an Nx-by-Nx interior grid, storing the result in Av.
void applyLaplacian(const std::vector<double>& v, std::vector<double>& Av, int Nx);

// Sum of the (up to 4) existing neighbor values of v at grid point (i, j),
// for use in point-wise updates such as Gauss-Seidel/SOR sweeps.
double applyLaplacianPoint(const std::vector<double>& v, int i, int j, int Nx);

// Dot product of two equal-length vectors
double dot(const std::vector<double>& a, const std::vector<double>& b);

// Euclidean (l2) norm of a vector
double norm2(const std::vector<double>& v);

#endif // LINALG_HPP
