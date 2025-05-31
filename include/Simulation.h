
#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
#include <memory> 
#include <cstdint> 





#include <SDL2/SDL_stdinc.h> 

#include "Vector3d.h"
#include "Camera.h"
#include "Particle.h"
#include "Octree.h"
#include "SPHGrid.h"
#include "SimConfig.h"

class Renderer; 
class BackgroundManager; 

class Simulation {
public:
    Simulation();
    ~Simulation();

    bool initialize(const char* window_title, int window_width, int window_height);
    void run();
    void cleanup();

private:
    std::unique_ptr<Renderer> renderer_ptr_;
    std::unique_ptr<BackgroundManager> background_manager_ptr_;

    Camera camera_;
    std::vector<Particle> particle_list_;
    Octree barnes_hut_tree_;
    SPHGrid sph_grid_object_;

    double simulation_initial_size_;
    double initial_camera_radius_;
    double initial_camera_azimuth_;
    double initial_camera_elevation_;
    double mouse_wheel_sensitivity_runtime_; 

    
    int num_dark_matter_particles_gui_;
    int num_gas_particles_gui_;
    float time_step_gui_;
    float G_gui_; 
    float H_SPH_gui_; 
    float K_EOS_gui_; 
    float EPSILON_SQUARED_gui_; 
    float MAX_EXPECTED_PARTICLE_SPEED_gui_; 
    float GAS_PARTICLE_JITTER_AMOUNT_gui_; 
    float theta_bh_gui_; 

    bool running_;
    bool simulation_paused_;
    double current_time_;
    Uint32 main_loop_counter_; 
    Uint32 total_physics_time_ms_; 
    int total_physics_steps_;
    Uint32 last_frame_ticks_; 

    bool mouse_left_button_down_;
    int mouse_prev_x_, mouse_prev_y_;

    void processEvents();
    void updateSimulationState(); 
    void stepPhysics(double dt);
    void renderFrame(float frame_delta_time_seconds);
    void setupGUI();
    void restartSimulation();
    void initializeParticlesInternal(int num_dm, int num_gas, double initial_size, double total_mass);

    void calculateForcesForChunk(
        size_t start_idx, size_t end_idx,
        const std::vector<Particle>& particles_read_only,
        std::vector<Vector3d>& forces_write,
        const Octree& tree, const SPHGrid& sph_grid_ref
    );
};

#endif 