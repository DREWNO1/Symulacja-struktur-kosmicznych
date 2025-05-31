#ifndef SIMCONFIG_H
#define SIMCONFIG_H




namespace SimConfig {
    
    const double PI = 3.14159265358979323846;

    
    inline double G = 1.0;

    
    inline double EPSILON_SQUARED = 0.01;

    
    inline double H_SPH = 10.0;

    
    inline double K_EOS = 0.1;

    
    inline double GAMMA_EOS = 5.0 / 3.0;

    
    inline float MAX_EXPECTED_PARTICLE_SPEED = 5.0f;

    
    inline float GAS_PARTICLE_JITTER_AMOUNT = 0.5f;

    
    inline double THETA_BARNES_HUT_SQUARED = 0.25;

    
    inline const double STARFIELD_DEPTH = 1500.0;
    inline const double DUST_FIELD_SCALE = 1.5; 
    inline const double DUST_MAX_SPEED = 0.05;

    inline const double MOUSE_SENSITIVITY = 0.005;
    
    
    
    inline const double MOUSE_WHEEL_SENSITIVITY_FACTOR = 1.0 / 20.0;
}

#endif 
