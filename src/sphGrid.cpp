#include "SPHGrid.h"

SPHGrid::SPHGrid()
    : min_bounds_(0,0,0),
      cell_size_(0.0),
      grid_dim_x_(0),
      grid_dim_y_(0),
      grid_dim_z_(0) {}

void SPHGrid::build(const std::vector<Particle>& particleList, double sph_smoothing_length) {
    if (particleList.empty() || sph_smoothing_length <= 0) {
        grid_dim_x_ = grid_dim_y_ = grid_dim_z_ = 0;
        cells_.clear();
        cell_size_ = 0.0;
        return;
    }

    Vector3d current_min_bounds(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Vector3d current_max_bounds(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
    bool gas_particles_exist = false;

    for (const auto& p : particleList) {
        if (p.type == ParticleType::GAS) {
            current_min_bounds.x = std::min(current_min_bounds.x, p.position.x);
            current_min_bounds.y = std::min(current_min_bounds.y, p.position.y);
            current_min_bounds.z = std::min(current_min_bounds.z, p.position.z);
            current_max_bounds.x = std::max(current_max_bounds.x, p.position.x);
            current_max_bounds.y = std::max(current_max_bounds.y, p.position.y);
            current_max_bounds.z = std::max(current_max_bounds.z, p.position.z);
            gas_particles_exist = true;
        }
    }

    if (!gas_particles_exist) {
        grid_dim_x_ = grid_dim_y_ = grid_dim_z_ = 0;
        cells_.clear();
        cell_size_ = 0.0;
        return;
    }

    cell_size_ = sph_smoothing_length; 
    double margin = cell_size_ * 2.0; 

    min_bounds_ = Vector3d(current_min_bounds.x - margin, current_min_bounds.y - margin, current_min_bounds.z - margin);
    Vector3d max_bounds_with_margin = Vector3d(current_max_bounds.x + margin, current_max_bounds.y + margin, current_max_bounds.z + margin);

    grid_dim_x_ = static_cast<int>(std::ceil((max_bounds_with_margin.x - min_bounds_.x) / cell_size_));
    grid_dim_y_ = static_cast<int>(std::ceil((max_bounds_with_margin.y - min_bounds_.y) / cell_size_));
    grid_dim_z_ = static_cast<int>(std::ceil((max_bounds_with_margin.z - min_bounds_.z) / cell_size_));

    grid_dim_x_ = std::max(1, grid_dim_x_);
    grid_dim_y_ = std::max(1, grid_dim_y_);
    grid_dim_z_ = std::max(1, grid_dim_z_);

    size_t total_cells = static_cast<size_t>(grid_dim_x_) * grid_dim_y_ * grid_dim_z_;
    if (cells_.size() != total_cells) {
        cells_.assign(total_cells, GridCell());
    } else {
        for (auto& cell : cells_) {
            cell.particle_indices.clear(); 
        }
    }

    for (size_t i = 0; i < particleList.size(); ++i) {
        const auto& p = particleList[i];
        if (p.type == ParticleType::GAS) {
            int cell_x, cell_y, cell_z;
            if (getCellCoordinates(p.position, cell_x, cell_y, cell_z)) {
                cells_[getFlatIndex(cell_x, cell_y, cell_z)].particle_indices.push_back(i);
            }
            
            
        }
    }
}

bool SPHGrid::getCellCoordinates(const Vector3d& position, int& out_cell_x, int& out_cell_y, int& out_cell_z) const {
    if (cell_size_ <= 0 || grid_dim_x_ == 0) return false; 

    out_cell_x = static_cast<int>(std::floor((position.x - min_bounds_.x) / cell_size_));
    out_cell_y = static_cast<int>(std::floor((position.y - min_bounds_.y) / cell_size_));
    out_cell_z = static_cast<int>(std::floor((position.z - min_bounds_.z) / cell_size_));

    
    if (out_cell_x >= 0 && out_cell_x < grid_dim_x_ &&
        out_cell_y >= 0 && out_cell_y < grid_dim_y_ &&
        out_cell_z >= 0 && out_cell_z < grid_dim_z_) {
        return true;
    }
    
    
    
    
    return false; 
}


size_t SPHGrid::getFlatIndex(int x, int y, int z) const {
    
    return static_cast<size_t>(x + y * grid_dim_x_ + z * grid_dim_x_ * grid_dim_y_);
}

const std::vector<size_t>& SPHGrid::getParticlesInCell(int cell_x, int cell_y, int cell_z) const {
    static const std::vector<size_t> empty_vector; 
    if (cell_x >= 0 && cell_x < grid_dim_x_ &&
        cell_y >= 0 && cell_y < grid_dim_y_ &&
        cell_z >= 0 && cell_z < grid_dim_z_ && isBuilt()) {
        return cells_[getFlatIndex(cell_x, cell_y, cell_z)].particle_indices;
    }
    return empty_vector;
}


std::vector<size_t> SPHGrid::getParticlesInNeighboringCells(int cell_x_base, int cell_y_base, int cell_z_base) const {
    std::vector<size_t> neighboring_particles;
    if (!isBuilt()) return neighboring_particles;

    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                int current_cell_x = cell_x_base + dx;
                int current_cell_y = cell_y_base + dy;
                int current_cell_z = cell_z_base + dz;

                
                if (current_cell_x >= 0 && current_cell_x < grid_dim_x_ &&
                    current_cell_y >= 0 && current_cell_y < grid_dim_y_ &&
                    current_cell_z >= 0 && current_cell_z < grid_dim_z_) {
                    
                    size_t flat_idx = getFlatIndex(current_cell_x, current_cell_y, current_cell_z);
                    if (flat_idx < cells_.size()) { 
                        const auto& particles_in_cell = cells_[flat_idx].particle_indices;
                        neighboring_particles.insert(neighboring_particles.end(), particles_in_cell.begin(), particles_in_cell.end());
                    }
                }
            }
        }
    }
    return neighboring_particles;
}

bool SPHGrid::isBuilt() const {
    return grid_dim_x_ > 0 && grid_dim_y_ > 0 && grid_dim_z_ > 0 && !cells_.empty();
}