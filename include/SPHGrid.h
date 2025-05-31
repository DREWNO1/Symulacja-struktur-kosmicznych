#ifndef SPHGRID_H
#define SPHGRID_H

#include <vector>
#include <cmath>     
#include <algorithm> 
#include <limits>    

#include "Vector3d.h"
#include "Particle.h"  
#include "SimConfig.h" 


struct GridCell {
    std::vector<size_t> particle_indices; 
};

class SPHGrid {
public:
    SPHGrid();

    
    void build(const std::vector<Particle>& particleList, double sph_smoothing_length);

    
    
    
    std::vector<size_t> getParticlesInNeighboringCells(int cell_x_base, int cell_y_base, int cell_z_base) const;
    
    
    const std::vector<size_t>& getParticlesInCell(int cell_x, int cell_y, int cell_z) const;


    
    
    
    bool getCellCoordinates(const Vector3d& position, int& out_cell_x, int& out_cell_y, int& out_cell_z) const;

    
    bool isBuilt() const;
    
    
    int getDimX() const { return grid_dim_x_; }
    int getDimY() const { return grid_dim_y_; }
    int getDimZ() const { return grid_dim_z_; }

private:
    std::vector<GridCell> cells_;         
    Vector3d min_bounds_;                 
    double cell_size_;                    
    int grid_dim_x_, grid_dim_y_, grid_dim_z_; 

    
    size_t getFlatIndex(int x, int y, int z) const;
};

#endif 