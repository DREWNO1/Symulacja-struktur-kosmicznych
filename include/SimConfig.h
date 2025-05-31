#ifndef SIMCONFIG_H
#define SIMCONFIG_H

// SimConfig.h
// Contains global simulation parameters and constants.

namespace SimConfig {
    // Mathematical constant PI
    const double PI = 3.14159265358979323846;

    // Gravitational constant (modifiable via ImGui in main.cpp)
    // This is declared as extern here because its value will be managed in main.cpp
    // and linked. If it were const, it could be defined directly here.
    // However, since it's modified at runtime, we'll define it in main.cpp
    // and just declare it here. For simplicity in a single-file-like structure
    // that was originally presented, we can define them here and let main.cpp modify them.
    // Let's keep them as global variables within the namespace.
    inline double G = 1.0;

    // Gravitational softening length squared, to prevent singularities.
    inline double EPSILON_SQUARED = 0.01;

    // SPH (Smoothed Particle Hydrodynamics) smoothing length.
    // Determines the radius of influence for SPH calculations.
    inline double H_SPH = 10.0;

    // Constant in the equation of state (EOS) for gas pressure.
    inline double K_EOS = 0.1;

    // Adiabatic index for gas, used in the equation of state.
    // Common value for a monatomic ideal gas.
    inline double GAMMA_EOS = 5.0 / 3.0;

    // Maximum expected particle speed.
    // Used for optimizing trail rendering dynamics. Adjust based on typical speeds in the simulation.
    inline float MAX_EXPECTED_PARTICLE_SPEED = 5.0f;

    // Max random displacement for rendering gas particles, adds a jitter effect.
    inline float GAS_PARTICLE_JITTER_AMOUNT = 0.5f;

    inline double THETA_BARNES_HUT_SQUARED = 0.25; 
}

#endif // SIMCONFIG_H
