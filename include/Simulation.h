// Simulation.h
#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>
#include <memory> // Dla std::unique_ptr
#include <cstdint> // Dla uint32_t

// Dołączenie nagłówka SDL dla typów SDL (jak Uint32, jeśli jest taka potrzeba)
// Alternatywnie, można użyć standardowych typów jak uint32_t z <cstdint>
// i rzutować w miejscach interakcji z SDL.
// Skoro błędy sugerują Uint32, dołączmy odpowiedni nagłówek SDL.
#include <SDL2/SDL_stdinc.h> // Zawiera definicję Uint32

#include "Vector3d.h"
#include "Camera.h"
#include "Particle.h"
#include "Octree.h"
#include "SPHGrid.h"
#include "SimConfig.h"

class Renderer; // Forward declaration
class BackgroundManager; // Forward declaration

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
    double mouse_wheel_sensitivity_runtime_; // Obliczane na podstawie simulationInitialSize

    // Pola dla wartości z GUI (zamiast globalnych statycznych lub w SimConfig)
    int num_dark_matter_particles_gui_;
    int num_gas_particles_gui_;
    float time_step_gui_;
    float G_gui_; // Odpowiada SimConfig::G
    float H_SPH_gui_; // Odpowiada SimConfig::H_SPH
    float K_EOS_gui_; // Odpowiada SimConfig::K_EOS
    float EPSILON_SQUARED_gui_; // Odpowiada SimConfig::EPSILON_SQUARED
    float MAX_EXPECTED_PARTICLE_SPEED_gui_; // Odpowiada SimConfig::MAX_EXPECTED_PARTICLE_SPEED
    float GAS_PARTICLE_JITTER_AMOUNT_gui_; // Odpowiada SimConfig::GAS_PARTICLE_JITTER_AMOUNT
    float theta_bh_gui_; // Przechowuje THETA (nie THETA^2) dla GUI

    bool running_;
    bool simulation_paused_;
    double current_time_;
    Uint32 main_loop_counter_; // Używamy typu z SDL
    Uint32 total_physics_time_ms_; // Używamy typu z SDL
    int total_physics_steps_;
    Uint32 last_frame_ticks_; // Używamy typu z SDL

    bool mouse_left_button_down_;
    int mouse_prev_x_, mouse_prev_y_;

    void processEvents();
    void updateSimulationState(); // Usunięto argument, bo deltaTime nie był tu potrzebny
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

#endif // SIMULATION_H