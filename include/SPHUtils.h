#ifndef SPHUTILS_H
#define SPHUTILS_H

#include "Vector3d.h"
#include "Particle.h"  
#include "SimConfig.h" 

namespace SPHUtils {

double kernelSPH(double r, double h);
double derivativeKernelSPH(double r, double h);
double calculatePressure(double density, double k_eos, double gamma_eos); 
Vector3d calculatePressureForceSymmetric(const Particle& p1, const Particle& p2, double h_sph, double k_eos, double gamma_eos); 

} 

#endif 