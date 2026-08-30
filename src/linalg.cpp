#include "linalg.hpp"
#include "mesh.hpp"
#include <cmath>

void applyLaplacian(const std::vector<double>& v, std::vector<double>& Av, int Nx) {
    for (int j = 0; j < Nx; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = index(i, j, Nx);
            double sum = 4.0 * v[idx];
            if (i > 0)      sum -= v[index(i - 1, j, Nx)];
            if (i < Nx - 1) sum -= v[index(i + 1, j, Nx)];
            if (j > 0)      sum -= v[index(i, j - 1, Nx)];
            if (j < Nx - 1) sum -= v[index(i, j + 1, Nx)];
            Av[idx] = sum;
        }
    }
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t k = 0; k < a.size(); ++k) sum += a[k] * b[k];
    return sum;
}

double norm2(const std::vector<double>& v) {
    return std::sqrt(dot(v, v));
}
