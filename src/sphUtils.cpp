#include "SPHUtils.h" 
#include <cmath>     




namespace SPHUtils {

    double kernelSPH(double r, double h) {
        if (r < 0.0 || r > h) {
            return 0.0;
        }
        if (h == 0.0) return 0.0;

        double h_sq = h * h;
        double r_sq = r * r;
        
        double h_ninth_calc = h_sq * h_sq * h_sq * h_sq * h; 
        if (h_ninth_calc == 0.0) return 0.0;

        double factor = 315.0 / (64.0 * SimConfig::PI * h_ninth_calc);
        double term = h_sq - r_sq;
        return factor * term * term * term; 
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
        return factor * term * term; 
    }

    
    double calculatePressure(double density, double k_eos, double gamma_eos) {
        if (density < 0.0) density = 0.0;
        return k_eos * std::pow(density, gamma_eos);
    }

    
    
    Vector3d calculatePressureForceSymmetric(const Particle& p1, const Particle& p2, double h_sph, double k_eos, double gamma_eos) {
        Vector3d r_ij = p2.position - p1.position;
        double dist = r_ij.length();

        if (dist >= 2.0 * h_sph || dist == 0.0) { 
            return Vector3d(0.0, 0.0, 0.0);
        }

        
        
        
        
        double pressure1 = p1.pressure;
        double pressure2 = p2.pressure;
        double rho1_sq = p1.density * p1.density;
        double rho2_sq = p2.density * p2.density;

        if (rho1_sq == 0.0 || rho2_sq == 0.0) {
            return Vector3d(0.0, 0.0, 0.0);
        }

        double pressure_term = (pressure1 / rho1_sq) + (pressure2 / rho2_sq);
        double grad_W = derivativeKernelSPH(dist, h_sph); 

        
        Vector3d unit_r_ij = r_ij.normalized(); 

        
        
        double scalar_factor = -p2.mass * pressure_term * grad_W;
        return unit_r_ij * scalar_factor;
    }

}