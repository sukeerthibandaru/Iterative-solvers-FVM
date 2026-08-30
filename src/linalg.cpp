#include "linalg.hpp"
#include "mesh.hpp"
#include <cmath>

void applyLaplacian(const std::vector<double>& v, std::vector<double>& Av, int Nx) {
    for (int j = 0; j < Nx; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = index(i, j, Nx);
            double sum = 4.0 * v[idx];
            if (i > 0)      sum -= v[index(i - 1, j, Nx)]; //left neighbor
            if (i < Nx - 1) sum -= v[index(i + 1, j, Nx)]; //right neighbor
            if (j > 0)      sum -= v[index(i, j - 1, Nx)]; //bottom neighbor 
            if (j < Nx - 1) sum -= v[index(i, j + 1, Nx)]; //top neighbor
            Av[idx] = sum;
        }
    }
}

double applyLaplacianPoint(const std::vector<double>& v, int i, int j, int Nx) {
    double sum = 0.0;
    if (i > 0)      sum += v[index(i - 1, j, Nx)]; //left neighbor
    if (i < Nx - 1) sum += v[index(i + 1, j, Nx)]; //right neighbor
    if (j > 0)      sum += v[index(i, j - 1, Nx)]; //bottom neighbor
    if (j < Nx - 1) sum += v[index(i, j + 1, Nx)]; //top neighbor
    return sum;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t k = 0; k < a.size(); ++k) sum += a[k] * b[k];
    return sum;
}

double norm2(const std::vector<double>& v) {
    return std::sqrt(dot(v, v));
}
