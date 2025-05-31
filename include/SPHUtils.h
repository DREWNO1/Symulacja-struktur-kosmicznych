#ifndef SPHUTILS_H
#define SPHUTILS_H

#include "Vector3d.h"
#include "Particle.h"  // Potrzebne dla calculatePressureForceSPHSymmetric
#include "SimConfig.h" // Dla SimConfig::PI, SimConfig::K_EOS, SimConfig::GAMMA_EOS, SimConfig::H_SPH

namespace SPHUtils {

double kernelSPH(double r, double h);
double derivativeKernelSPH(double r, double h);
double calculatePressure(double density, double k_eos, double gamma_eos); // Zamiast calculatePressureSPH, żeby było bardziej generyczne
Vector3d calculatePressureForceSymmetric(const Particle& p1, const Particle& p2, double h_sph, double k_eos, double gamma_eos); // Przemianowane i bardziej generyczne

} // namespace SPHUtils

#endif // SPHUTILS_H