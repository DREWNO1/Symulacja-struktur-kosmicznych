#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vector3d.h" 
#include <deque>      
#include <cstddef>    



enum ParticleType {
    DARK_MATTER,
    GAS
    
};






struct Particle {
    Vector3d position;
    Vector3d velocity;
    double mass;
    ParticleType type;
    double density;     
    double pressure;    
    std::deque<Vector3d> trail_history; 
    size_t id;          

    
    Particle();

    
    Particle(size_t p_id, 
             ParticleType p_type = ParticleType::DARK_MATTER, 
             double p_mass = 0.0,
             const Vector3d& p_pos = Vector3d(0,0,0),
             const Vector3d& p_vel = Vector3d(0,0,0));

    
    
    
};

#endif 