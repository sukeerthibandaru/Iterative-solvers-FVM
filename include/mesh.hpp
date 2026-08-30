#ifndef MESH_HPP
#define MESH_HPP

#include <vector>

struct Mesh {
    int N;                      // Total grid points per dimension (including boundaries)
    int Nx;                     // Interior grid points per dimension (N - 2)
    double h;                   // Uniform grid spacing h = 1 / (N - 1)
    std::vector<double> f;      // Discretized Right-Hand Side (RHS) source vector
};

// Helper function to map 2D interior cell indices (i, j) to 1D flat vector index
inline int index(int i, int j, int Nx) {
    return j * Nx + i;
}

// Function prototype to build and initialize the 2D grid and source vector
Mesh createMesh(int N);

#endif // MESH_HPP