#include "Particle.h"

Particle::Particle() 
    : position(0.0, 0.0, 0.0),
      velocity(0.0, 0.0, 0.0),
      mass(0.0),
      type(ParticleType::DARK_MATTER),
      density(0.0),
      pressure(0.0),
      id(0)
{}

Particle::Particle(size_t p_id, 
                   ParticleType p_type, 
                   double p_mass,
                   const Vector3d& p_pos,
                   const Vector3d& p_vel)
    : position(p_pos),
      velocity(p_vel),
      mass(p_mass),
      type(p_type),
      density(0.0),
      pressure(0.0),
      id(p_id)
{}