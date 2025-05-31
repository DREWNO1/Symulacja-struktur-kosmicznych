#ifndef OCTREE_H
#define OCTREE_H

#include "Vector3d.h"
#include "SimConfig.h" // Potrzebne dla SimConfig::THETA_BARNES_HUT_SQUARED i SimConfig::EPSILON_SQUARED
#include <vector>
#include <memory> // Dla std::unique_ptr
#include <limits> // Dla std::numeric_limits
#include <algorithm> // Dla std::min, std::max

// --- Forward declaration for Particle struct ---
// Zakładamy, że Particle będzie zdefiniowane gdzie indziej (np. w pliku z główną symulacją lub particle.h)
struct Particle; // Pełna definicja Particle musi być dostępna dla octree.cpp

// --- Barnes-Hut Octree Structures ---
const int MAX_PARTICLES_PER_OCTREE_NODE = 1; // Zwykle 1 dla klasycznego Barnes-Hut

struct OctreeNode {
    Vector3d center;        // Środek geometryczny węzła
    double half_width;      // Połowa szerokości węzła (zakładamy sześciany)
    std::unique_ptr<OctreeNode> children[8]; // Dzieci węzła (oktanty)
    
    // Dane o masie
    double total_mass;
    Vector3d center_of_mass;
    
    // Cząstki w tym węźle
    std::vector<size_t> particle_indices_in_node; // Przechowuje indeksy cząstek *tylko* jeśli jest liściem i ma > 0 cząstek
                                                 // lub tymczasowo podczas budowy, zanim się podzieli.
    bool is_leaf;
    int num_particles_in_subtree; // Liczba cząstek w tym węźle i jego poddrzewach

    OctreeNode(const Vector3d& c, double hw);

    // Zwraca, do którego oktantu (0-7) należy dana pozycja względem środka tego węzła
    int getOctantContainingPoint(const Vector3d& point) const;

    // Dzieli węzeł na 8 dzieci
    void subdivide();
};

class Octree {
public:
    std::unique_ptr<OctreeNode> root;
    const std::vector<Particle>* all_particles_ptr; // Wskaźnik do globalnej listy cząstek

    Octree();

    void build(const std::vector<Particle>& particleList);
    void insertParticle(OctreeNode* node, size_t particle_idx);
    void computeMassDistribution(OctreeNode* node);
    Vector3d calculateForce(const Particle& target_particle, size_t target_particle_id) const;

private:
    Vector3d calculateForceRecursive(const Particle& target_particle, size_t target_particle_id, const OctreeNode* node) const;
    // Funkcja calculateGravityForce jest globalna lub statyczna, więc nie jest częścią klasy Octree,
    // ale jest używana przez Octree. Jej deklaracja/definicja musi być dostępna.
    // Jeśli chcemy ją przenieść jako metodę prywatną (np. statyczną), to tutaj byłaby deklaracja.
    // Na razie zakładamy, że jest dostępna globalnie, jak w oryginalnym main.cpp.
};

#endif // OCTREE_H