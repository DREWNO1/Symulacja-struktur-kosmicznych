#ifndef SIMCONFIG_H
#define SIMCONFIG_H

// SimConfig.h
// Contains global simulation parameters and constants.

namespace SimConfig {
    // Mathematical constant PI
    const double PI = 3.14159265358979323846;

    // Gravitational constant
    inline double G = 1.0;

    // Gravitational softening length squared
    inline double EPSILON_SQUARED = 0.01;

    // SPH (Smoothed Particle Hydrodynamics) smoothing length
    inline double H_SPH = 10.0;

    // Constant in the equation of state (EOS) for gas pressure
    inline double K_EOS = 0.1;

    // Adiabatic index for gas
    inline double GAMMA_EOS = 5.0 / 3.0;

    // Maximum expected particle speed (for trail rendering dynamics)
    inline float MAX_EXPECTED_PARTICLE_SPEED = 5.0f;

    // Max random displacement for rendering gas particles
    inline float GAS_PARTICLE_JITTER_AMOUNT = 0.5f;

    // Barnes-Hut opening angle parameter squared (theta^2)
    inline double THETA_BARNES_HUT_SQUARED = 0.25;

    // --- NOWE STAŁE PRZENIESIONE Z MAINA ---
    inline const double STARFIELD_DEPTH = 1500.0;
    inline const double DUST_FIELD_SCALE = 1.5; // Było używane w main.cpp z simulationInitialSize
    inline const double DUST_MAX_SPEED = 0.05;

    inline const double MOUSE_SENSITIVITY = 0.005;
    // MOUSE_WHEEL_SENSITIVITY zależało od simulationInitialSize,
    // więc lepiej zostawić je jako obliczane w Simulation::initialize lub przekazywać.
    // Lub zdefiniować tu jako współczynnik:
    inline const double MOUSE_WHEEL_SENSITIVITY_FACTOR = 1.0 / 20.0;
}

#endif // SIMCONFIG_H
