#include "octree.h"
#include "Particle.h" // Musi zawierać definicję struktury Particle
                      // (zakładam, że stworzysz taki plik lub Particle jest w innym dostępnym nagłówku)

// --- Implementacja OctreeNode ---
OctreeNode::OctreeNode(const Vector3d& c, double hw)
    : center(c), half_width(hw), total_mass(0.0), center_of_mass(0,0,0), is_leaf(true), num_particles_in_subtree(0) {
    for (int i = 0; i < 8; ++i) {
        children[i] = nullptr;
    }
}

int OctreeNode::getOctantContainingPoint(const Vector3d& point) const {
    int octant = 0;
    if (point.x >= center.x) octant |= 1; // Prawo
    if (point.y >= center.y) octant |= 2; // Góra
    if (point.z >= center.z) octant |= 4; // Przód
    return octant;
}

void OctreeNode::subdivide() {
    double child_half_width = half_width / 2.0;
    for (int i = 0; i < 8; ++i) {
        Vector3d child_center = center;
        child_center.x += (i & 1) ? child_half_width : -child_half_width;
        child_center.y += (i & 2) ? child_half_width : -child_half_width;
        child_center.z += (i & 4) ? child_half_width : -child_half_width;
        children[i] = std::make_unique<OctreeNode>(child_center, child_half_width);
    }
    is_leaf = false; // Już nie jest liściem
}

// --- Implementacja Octree ---
Octree::Octree() : root(nullptr), all_particles_ptr(nullptr) {}

void Octree::build(const std::vector<Particle>& particleList) {
    all_particles_ptr = &particleList;
    if (particleList.empty()) {
        root = nullptr;
        return;
    }

    // 1. Znajdź granice dla korzenia drzewa
    Vector3d min_coord(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Vector3d max_coord(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());

    for (const auto& p : particleList) {
        min_coord.x = std::min(min_coord.x, p.position.x);
        min_coord.y = std::min(min_coord.y, p.position.y);
        min_coord.z = std::min(min_coord.z, p.position.z);
        max_coord.x = std::max(max_coord.x, p.position.x);
        max_coord.y = std::max(max_coord.y, p.position.y);
        max_coord.z = std::max(max_coord.z, p.position.z);
    }

    Vector3d center = (min_coord + max_coord) * 0.5;
    double half_width = 0.5 * std::max({max_coord.x - min_coord.x, max_coord.y - min_coord.y, max_coord.z - min_coord.z});
    half_width = std::max(half_width, 1e-5); // Unikaj zerowej szerokości

    root = std::make_unique<OctreeNode>(center, half_width);

    // 2. Wstaw wszystkie cząstki
    for (size_t i = 0; i < particleList.size(); ++i) {
        insertParticle(root.get(), i);
    }

    // 3. Oblicz rozkład masy
    computeMassDistribution(root.get());
}

void Octree::insertParticle(OctreeNode* node, size_t particle_idx) {
    if (!node) return;

    node->num_particles_in_subtree++; 

    if (node->is_leaf) {
        if (node->particle_indices_in_node.empty()) { 
            node->particle_indices_in_node.push_back(particle_idx);
        } else {
            size_t existing_particle_idx = node->particle_indices_in_node[0];
            node->particle_indices_in_node.clear(); 
            
            node->subdivide(); 

            int octant_existing = node->getOctantContainingPoint((*all_particles_ptr)[existing_particle_idx].position);
            insertParticle(node->children[octant_existing].get(), existing_particle_idx);

            int octant_new = node->getOctantContainingPoint((*all_particles_ptr)[particle_idx].position);
            insertParticle(node->children[octant_new].get(), particle_idx);
        }
    } else { 
        int octant = node->getOctantContainingPoint((*all_particles_ptr)[particle_idx].position);
        insertParticle(node->children[octant].get(), particle_idx);
    }
}
    
void Octree::computeMassDistribution(OctreeNode* node) {
    if (!node) return;

    if (node->is_leaf) {
        if (!node->particle_indices_in_node.empty()) {
            for(size_t p_idx : node->particle_indices_in_node) {
                 node->total_mass += (*all_particles_ptr)[p_idx].mass;
                 node->center_of_mass = node->center_of_mass + (*all_particles_ptr)[p_idx].position * (*all_particles_ptr)[p_idx].mass;
            }
            if (node->total_mass > 0) {
                node->center_of_mass = node->center_of_mass / node->total_mass;
            } else {
                 node->center_of_mass = node->center; 
            }
        } else { 
            node->total_mass = 0.0;
            node->center_of_mass = node->center;
        }
    } else { 
        node->total_mass = 0.0;
        node->center_of_mass = Vector3d(0,0,0);
        for (int i = 0; i < 8; ++i) {
            if (node->children[i]) {
                computeMassDistribution(node->children[i].get());
                node->total_mass += node->children[i]->total_mass;
                node->center_of_mass = node->center_of_mass + node->children[i]->center_of_mass * node->children[i]->total_mass;
            }
        }
        if (node->total_mass > 0) {
            node->center_of_mass = node->center_of_mass / node->total_mass;
        } else {
            node->center_of_mass = node->center; 
        }
    }
}

Vector3d Octree::calculateForce(const Particle& target_particle, size_t target_particle_id) const {
    if (!root) return Vector3d(0,0,0);
    return calculateForceRecursive(target_particle, target_particle_id, root.get());
}

// Funkcja pomocnicza do obliczania siły grawitacji między dwiema cząstkami.
// Została przeniesiona tutaj jako funkcja statyczna lub wolna funkcja w tym pliku .cpp,
// aby uniknąć zależności od globalnej funkcji w main.cpp, jeśli jest to pożądane.
// Jeśli jest to globalna funkcja, ten blok nie jest potrzebny, a jedynie jej deklaracja
// (lub upewnienie się, że odpowiedni nagłówek jest dołączony).
// Zakładając, że definicja `Particle` zawiera `position` i `mass`.
// Ta funkcja powinna być zdefiniowana przed jej użyciem w `calculateForceRecursive`
// lub zadeklarowana w `octree.h` (jeśli ma być dostępna również na zewnątrz).
// Dla uproszczenia, kopiuję ją tutaj jako funkcję statyczną wewnątrz pliku.
// Można ją też umieścić w przestrzeni nazw anonimowej.
namespace { // Anonimowa przestrzeń nazw ogranicza widoczność do tego pliku
    Vector3d calculateDirectGravityForce(const Particle& p1, const Particle& p2) {
        Vector3d r_vec = p2.position - p1.position;
        double dist_sq = r_vec.lengthSquared();

        if (dist_sq == 0.0 && SimConfig::EPSILON_SQUARED == 0.0) return Vector3d(0.0, 0.0, 0.0);

        double dist_sq_softened = dist_sq + SimConfig::EPSILON_SQUARED;
        if (dist_sq_softened == 0.0) return Vector3d(0.0, 0.0, 0.0);

        double force_mag = SimConfig::G * p1.mass * p2.mass / dist_sq_softened;
        
        // Aby uniknąć normalizacji wektora zerowego, jeśli dist_sq == 0 ale EPSILON_SQUARED > 0
        if (dist_sq == 0.0) return Vector3d(0.0,0.0,0.0); // Siła powinna być zerowa w tym przypadku symetrii

        return r_vec.normalized() * force_mag;
    }
} // namespace


Vector3d Octree::calculateForceRecursive(const Particle& target_particle, size_t target_particle_id, const OctreeNode* node) const {
    if (!node || node->num_particles_in_subtree == 0) {
        return Vector3d(0,0,0);
    }

    if (node->is_leaf) {
        Vector3d force(0,0,0);
        for(size_t p_idx_in_leaf : node->particle_indices_in_node) {
            // Używamy target_particle.id (zakładając, że Particle ma pole 'id')
            // lub porównujemy wskaźniki/indeksy, jeśli ID nie ma.
            // Oryginalny kod używał target_particle_id, który był indeksem.
            // Zakładamy, że all_particles_ptr[p_idx_in_leaf].id istnieje.
            if (p_idx_in_leaf != target_particle_id) { 
                 force = force + calculateDirectGravityForce(target_particle, (*all_particles_ptr)[p_idx_in_leaf]);
            }
        }
        return force;
    }

    double s = node->half_width * 2.0;
    double d_squared = (target_particle.position - node->center_of_mass).lengthSquared();

    if (d_squared == 0.0) {
        // Cząstka jest w środku masy węzła. Musimy zejść głębiej, aby uniknąć dzielenia przez zero
        // i uzyskać dokładniejszy wynik, chyba że EPSILON_SQUARED jest wystarczająco duże.
        // W tej sytuacji s/d jest nieskończone, więc warunek theta nie zostanie spełniony.
        // Powinno się zawsze rekurencyjnie schodzić głębiej.
        Vector3d force_sum(0,0,0);
        for(int i=0; i<8; ++i) {
            if(node->children[i] && node->children[i]->num_particles_in_subtree > 0) {
                force_sum = force_sum + calculateForceRecursive(target_particle, target_particle_id, node->children[i].get());
            }
        }
        return force_sum;
    }

    if ((s * s / d_squared) < SimConfig::THETA_BARNES_HUT_SQUARED) {
        if (node->total_mass == 0.0) return Vector3d(0,0,0);

        Vector3d direction = node->center_of_mass - target_particle.position;
        // d_squared jest już obliczone. Softening jest dodawany tutaj.
        double dist_sq_with_softening = d_squared + SimConfig::EPSILON_SQUARED;
        
        if (dist_sq_with_softening == 0) return Vector3d(0,0,0); 

        double force_magnitude = SimConfig::G * target_particle.mass * node->total_mass / dist_sq_with_softening;
        Vector3d unit_direction = direction.normalized(); 
        return unit_direction * force_magnitude;

    } else {
        Vector3d force_sum(0,0,0);
        for (int i = 0; i < 8; ++i) {
            if (node->children[i] && node->children[i]->num_particles_in_subtree > 0) {
                force_sum = force_sum + calculateForceRecursive(target_particle, target_particle_id, node->children[i].get());
            }
        }
        return force_sum;
    }
}