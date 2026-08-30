#include "mesh.hpp"
#include <cmath>

const double PI = std::acos(-1.0);

Mesh createMesh(int N) {
    Mesh mesh;
    mesh.N = N;
    mesh.Nx = N - 2;
    mesh.h = 1.0 / (N - 1);
    
    int size = mesh.Nx * mesh.Nx;
    mesh.f.resize(size, 0.0);

    // Populate Right-Hand Side (RHS) source term f(x, y) * h^2
    for (int j = 0; j < mesh.Nx; ++j) {
        for (int i = 0; i < mesh.Nx; ++i) {
            double x = (i + 1) * mesh.h;
            double y = (j + 1) * mesh.h;

            mesh.f[index(i, j, mesh.Nx)] = mesh.h * mesh.h * (
                -2.0 * PI * PI * std::cos(2.0 * PI * x) * std::pow(std::sin(PI * y), 2)
                -2.0 * PI * PI * std::pow(std::sin(PI * x), 2) * std::cos(2.0 * PI * y)  
            );
        }
    }

    return mesh;
}