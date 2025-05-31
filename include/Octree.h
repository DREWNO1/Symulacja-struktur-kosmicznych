#ifndef OCTREE_H
#define OCTREE_H

#include "Vector3d.h"
#include "SimConfig.h" 
#include <vector>
#include <memory> 
#include <limits> 
#include <algorithm> 



struct Particle; 


const int MAX_PARTICLES_PER_OCTREE_NODE = 1; 

struct OctreeNode {
    Vector3d center;        
    double half_width;      
    std::unique_ptr<OctreeNode> children[8]; 
    
    
    double total_mass;
    Vector3d center_of_mass;
    
    
    std::vector<size_t> particle_indices_in_node; 
                                                 
    bool is_leaf;
    int num_particles_in_subtree; 

    OctreeNode(const Vector3d& c, double hw);

    
    int getOctantContainingPoint(const Vector3d& point) const;

    
    void subdivide();
};

class Octree {
public:
    std::unique_ptr<OctreeNode> root;
    const std::vector<Particle>* all_particles_ptr; 

    Octree();

    void build(const std::vector<Particle>& particleList);
    void insertParticle(OctreeNode* node, size_t particle_idx);
    void computeMassDistribution(OctreeNode* node);
    Vector3d calculateForce(const Particle& target_particle, size_t target_particle_id) const;

private:
    Vector3d calculateForceRecursive(const Particle& target_particle, size_t target_particle_id, const OctreeNode* node) const;
    
    
    
    
};

#endif 