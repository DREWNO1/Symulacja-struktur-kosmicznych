#include "SPHUtils.h" // Zakładamy, że SPHUtils.h jest już zdefiniowany
#include <cmath>     // Dla std::pow

// Definicje z main.cpp
// double SimConfig::PI = 3.14159265358979323846; // PI jest teraz w SimConfig.h

namespace SPHUtils {

    double kernelSPH(double r, double h) {
        if (r < 0.0 || r > h) {
            return 0.0;
        }
        if (h == 0.0) return 0.0;

        double h_sq = h * h;
        double r_sq = r * r;
        // double h_ninth = std::pow(h, 9.0); // Można zoptymalizować, aby uniknąć pow()
        double h_ninth_calc = h_sq * h_sq * h_sq * h_sq * h; // h^9
        if (h_ninth_calc == 0.0) return 0.0;

        double factor = 315.0 / (64.0 * SimConfig::PI * h_ninth_calc);
        double term = h_sq - r_sq;
        return factor * term * term * term; // (h^2 - r^2)^3
    }

    double derivativeKernelSPH(double r, double h) {
        if (r < 0.0 || r > h) {
            return 0.0;
        }
        if (h == 0.0) return 0.0;

        double h_sixth = std::pow(h, 6.0);
        if (h_sixth == 0.0) return 0.0;

        double factor = -45.0 / (SimConfig::PI * h_sixth);
        double term = h - r;
        return factor * term * term; // factor * (h-r)^2
    }

    // Zmieniona nazwa i sygnatura w stosunku do main.cpp
    double calculatePressure(double density, double k_eos, double gamma_eos) {
        if (density < 0.0) density = 0.0;
        return k_eos * std::pow(density, gamma_eos);
    }

    // Zmieniona nazwa i sygnatura w stosunku do main.cpp
    // Przyjmuje parametry konfiguracyjne bezpośrednio.
    Vector3d calculatePressureForceSymmetric(const Particle& p1, const Particle& p2, double h_sph, double k_eos, double gamma_eos) {
        Vector3d r_ij = p2.position - p1.position;
        double dist = r_ij.length();

        if (dist >= 2.0 * h_sph || dist == 0.0) { // Zgodnie z oryginalnym kodem, 2.0*h_sph
            return Vector3d(0.0, 0.0, 0.0);
        }

        // Ciśnienia są polami cząstek, więc je używamy. Gęstości również.
        // Jeśli ciśnienia nie są jeszcze obliczone, ta funkcja nie powinna być wołana
        // lub powinna przyjmować gęstości i obliczać ciśnienia w locie.
        // Zakładam, że p1.pressure i p2.pressure są już obliczone.
        double pressure1 = p1.pressure;
        double pressure2 = p2.pressure;
        double rho1_sq = p1.density * p1.density;
        double rho2_sq = p2.density * p2.density;

        if (rho1_sq == 0.0 || rho2_sq == 0.0) {
            return Vector3d(0.0, 0.0, 0.0);
        }

        double pressure_term = (pressure1 / rho1_sq) + (pressure2 / rho2_sq);
        double grad_W = derivativeKernelSPH(dist, h_sph); // dW/dr

        // Wektor jednostkowy r_ij / dist (normalizacja r_ij)
        Vector3d unit_r_ij = r_ij.normalized(); // Zabezpieczone przed dzieleniem przez zero w Vector3d::normalized

        // F_i = -m_j * sum_j (P_i/rho_i^2 + P_j/rho_j^2) * grad_W_ij * unit_r_ij
        // Siła działająca na p1 od p2
        double scalar_factor = -p2.mass * pressure_term * grad_W;
        return unit_r_ij * scalar_factor;
    }

}