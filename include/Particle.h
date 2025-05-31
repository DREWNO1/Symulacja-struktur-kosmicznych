#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vector3d.h" // Zakładamy, że Vector3d.h istnieje i jest dostępny
#include <deque>      // Dla std::deque (używane dla historii śladu)
#include <cstddef>    // Dla size_t

// Definicja ParticleType, która była w main.cpp, powinna zostać przeniesiona tutaj
// lub do wspólnego pliku nagłówkowego z typami, jeśli Particle ma z niej korzystać.
enum ParticleType {
    DARK_MATTER,
    GAS
    // Można dodać UNKNOWN lub DEFAULT, jeśli to potrzebne
};

// Maksymalna historia śladu, jeśli chcemy to uczynić konfigurowalne per cząstkę
// lub przenieść do SimConfig. W oryginalnym kodzie było to globalne.
// const int MAX_TRAIL_HISTORY_PARTICLE = 10;


struct Particle {
    Vector3d position;
    Vector3d velocity;
    double mass;
    ParticleType type;
    double density;     // Zwykle obliczane dla cząstek gazu w SPH
    double pressure;    // Zwykle obliczane dla cząstek gazu w SPH
    std::deque<Vector3d> trail_history; // Historia pozycji dla rysowania śladów
    size_t id;          // Unikalne ID dla każdej cząstki

    // Deklaracja konstruktora domyślnego
    Particle();

    // Przykładowy konstruktor z parametrami (można go rozbudować)
    Particle(size_t p_id, 
             ParticleType p_type = ParticleType::DARK_MATTER, 
             double p_mass = 0.0,
             const Vector3d& p_pos = Vector3d(0,0,0),
             const Vector3d& p_vel = Vector3d(0,0,0));

    // Można tu dodać deklaracje innych metod, np.:
    // void updateTrail(const Vector3d& new_position, size_t max_length);
    // void resetState();
};

#endif // PARTICLE_H