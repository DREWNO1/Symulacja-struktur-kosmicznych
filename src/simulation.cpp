// Simulation.cpp
#include "Simulation.h"
#include "RandomUtils.h"
#include "SPHUtils.h"
#include "Renderer.h" // Dołączony dla std::unique_ptr<Renderer>
#include "BackgroundManager.h" // Dołączony dla std::unique_ptr<BackgroundManager>
#include <thread>      // <<< BRAKUJĄCY INCLUDE



// Dołączenia dla ImGui (jeśli logika ImGui jest w Simulation)
#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_sdlrenderer2.h"

// Stałe, które były globalne w main.cpp, a nie pasują do SimConfig
const int MAX_TRAIL_HISTORY_SIM = 10; // Przemianowane, aby uniknąć konfliktu z ewentualną stałą w Particle.h

Simulation::Simulation()
    : renderer_ptr_(nullptr), background_manager_ptr_(nullptr),
      simulation_initial_size_(200.0),
      initial_camera_radius_(simulation_initial_size_ * 2.25),
      initial_camera_azimuth_(0.0),
      initial_camera_elevation_(0.0),
      mouse_wheel_sensitivity_runtime_(0.0), // Obliczone w initialize
      // Inicjalizacja pól _gui_ wartościami z SimConfig
      num_dark_matter_particles_gui_(500), // Domyślne wartości, mogą być nadpisane
      num_gas_particles_gui_(2000),
      time_step_gui_(0.005f),
      G_gui_(static_cast<float>(SimConfig::G)),
      H_SPH_gui_(static_cast<float>(SimConfig::H_SPH)),
      K_EOS_gui_(static_cast<float>(SimConfig::K_EOS)),
      EPSILON_SQUARED_gui_(static_cast<float>(SimConfig::EPSILON_SQUARED)),
      MAX_EXPECTED_PARTICLE_SPEED_gui_(SimConfig::MAX_EXPECTED_PARTICLE_SPEED),
      GAS_PARTICLE_JITTER_AMOUNT_gui_(SimConfig::GAS_PARTICLE_JITTER_AMOUNT),
      theta_bh_gui_(static_cast<float>(std::sqrt(SimConfig::THETA_BARNES_HUT_SQUARED))),
      running_(false), simulation_paused_(false), current_time_(0.0),
      main_loop_counter_(0), total_physics_time_ms_(0), total_physics_steps_(0), // Inicjalizacja pól
      last_frame_ticks_(0), mouse_left_button_down_(false), mouse_prev_x_(0), mouse_prev_y_(0)
{
}

Simulation::~Simulation() {
    cleanup();
}


bool Simulation::initialize(const char* window_title, int window_width, int window_height) {
    renderer_ptr_ = std::make_unique<Renderer>();
    if (!renderer_ptr_->initialize(window_title, window_width, window_height)) {
        return false;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplSDL2_InitForSDLRenderer(renderer_ptr_->getSDLWindow(), renderer_ptr_->getSDLRenderer());
    ImGui_ImplSDLRenderer2_Init(renderer_ptr_->getSDLRenderer());

    camera_.initializeSpherical(initial_camera_radius_, initial_camera_azimuth_, initial_camera_elevation_);
    camera_.lookAt = Vector3d(0, 0, 0);
    camera_.up = Vector3d(0, 1, 0);
    camera_.fovY_degrees = 60.0;
    camera_.nearClip = simulation_initial_size_ / 100.0;
    camera_.farClip = simulation_initial_size_ * 10.0;
    camera_.lerp_speed = 0.08f;

    mouse_wheel_sensitivity_runtime_ = simulation_initial_size_ * SimConfig::MOUSE_WHEEL_SENSITIVITY_FACTOR;

    background_manager_ptr_ = std::make_unique<BackgroundManager>();
    background_manager_ptr_->initialize(simulation_initial_size_, 10.0, SimConfig::DUST_FIELD_SCALE);

    initializeParticlesInternal(num_dark_matter_particles_gui_, num_gas_particles_gui_, simulation_initial_size_, 1e5);

    current_time_ = 0.0;
    running_ = true;
    simulation_paused_ = false;
    main_loop_counter_ = 0; // Poprawne użycie pola klasy
    total_physics_time_ms_ = 0; // Poprawne użycie pola klasy
    total_physics_steps_ = 0;
    last_frame_ticks_ = SDL_GetTicks(); // Poprawne użycie pola klasy

    return true;
}

void Simulation::cleanup() {
    ImGui_ImplSDLRenderer2_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    if (ImGui::GetCurrentContext()) { // Sprawdź, czy kontekst istnieje przed zniszczeniem
        ImGui::DestroyContext();
    }

    if (renderer_ptr_) {
        renderer_ptr_->shutdown();
    }
}

void Simulation::run() {
    if (!running_) return;

    while (running_) {
        Uint32 current_ticks = SDL_GetTicks();
        float frame_delta_time_seconds = (current_ticks - last_frame_ticks_) / 1000.0f;
        if (frame_delta_time_seconds <= 0) frame_delta_time_seconds = 1.0f / 60.0f;
        if (frame_delta_time_seconds > 0.1f) frame_delta_time_seconds = 0.1f;
        last_frame_ticks_ = current_ticks;

        processEvents();
        updateSimulationState(); 
        
        if (!simulation_paused_) {
            if (background_manager_ptr_) {
                 background_manager_ptr_->update(frame_delta_time_seconds, simulation_initial_size_, SimConfig::DUST_FIELD_SCALE);
            }
            if (time_step_gui_ > 0) {
                Uint32 sim_step_start_ticks = SDL_GetTicks();
                stepPhysics(static_cast<double>(time_step_gui_));
                Uint32 sim_step_duration_ms = SDL_GetTicks() - sim_step_start_ticks;
                total_physics_time_ms_ += sim_step_duration_ms;
                total_physics_steps_++;
                current_time_ += static_cast<double>(time_step_gui_);
            }
        }
        main_loop_counter_++;

        renderFrame(frame_delta_time_seconds);

        Uint32 frame_render_duration_ms = SDL_GetTicks() - current_ticks;
        const Uint32 TARGET_FRAME_DURATION_MS = 1000 / 60;
        if (frame_render_duration_ms < TARGET_FRAME_DURATION_MS) {
            SDL_Delay(TARGET_FRAME_DURATION_MS - frame_render_duration_ms);
        }
    }
}

void Simulation::processEvents() {
    SDL_Event event;
    ImGuiIO& io = ImGui::GetIO();

    while (SDL_PollEvent(&event)) {
        ImGui_ImplSDL2_ProcessEvent(&event);

        bool imgui_captures_mouse = io.WantCaptureMouse;
        bool imgui_captures_keyboard = io.WantCaptureKeyboard;

        if (event.type == SDL_QUIT) {
            running_ = false;
        }
        // Window resize jest obsługiwane przez Renderer::beginFrame()

        if (!imgui_captures_mouse) {
            if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_LEFT) {
                mouse_left_button_down_ = true;
                SDL_GetMouseState(&mouse_prev_x_, &mouse_prev_y_);
            }
            if (event.type == SDL_MOUSEBUTTONUP && event.button.button == SDL_BUTTON_LEFT) {
                mouse_left_button_down_ = false;
            }
            if (event.type == SDL_MOUSEMOTION && mouse_left_button_down_) {
                int current_mouse_x, current_mouse_y;
                SDL_GetMouseState(&current_mouse_x, &current_mouse_y);
                int delta_x = current_mouse_x - mouse_prev_x_;
                int delta_y = current_mouse_y - mouse_prev_y_;
                camera_.target_azimuth -= static_cast<double>(delta_x) * SimConfig::MOUSE_SENSITIVITY;
                camera_.target_elevation -= static_cast<double>(delta_y) * SimConfig::MOUSE_SENSITIVITY;
                mouse_prev_x_ = current_mouse_x;
                mouse_prev_y_ = current_mouse_y;
            }
            if (event.type == SDL_MOUSEWHEEL) {
                if (event.wheel.y > 0) camera_.target_radius -= mouse_wheel_sensitivity_runtime_;
                else if (event.wheel.y < 0) camera_.target_radius += mouse_wheel_sensitivity_runtime_;
            }
        }

        if (!imgui_captures_keyboard) {
            if (event.type == SDL_KEYDOWN) {
                switch (event.key.keysym.sym) {
                    case SDLK_r: restartSimulation(); break;
                    case SDLK_p: simulation_paused_ = !simulation_paused_; break;
                    case SDLK_LEFT: camera_.target_azimuth -= camera_.lerp_speed * 5.0f; break;
                    case SDLK_RIGHT: camera_.target_azimuth += camera_.lerp_speed * 5.0f; break;
                    case SDLK_UP: camera_.target_elevation += camera_.lerp_speed * 5.0f; break;
                    case SDLK_DOWN: camera_.target_elevation -= camera_.lerp_speed * 5.0f; break;
                    case SDLK_a: camera_.target_radius -= mouse_wheel_sensitivity_runtime_; break;
                    case SDLK_z: camera_.target_radius += mouse_wheel_sensitivity_runtime_; break;
                    default: break;
                }
            }
        }
        camera_.target_elevation = std::max(-SimConfig::PI / 2.0 + 0.02, std::min(SimConfig::PI / 2.0 - 0.02, camera_.target_elevation));
        camera_.target_radius = std::max(camera_.nearClip * 3.0, camera_.target_radius);
    }
}


void Simulation::updateSimulationState() { // Usunięto frame_delta_time_seconds
    camera_.smoothUpdate();

    if (!simulation_paused_ && time_step_gui_ > 0) {
        SimConfig::G = G_gui_; // Użyj pól klasy Simulation
        SimConfig::H_SPH = H_SPH_gui_;
        SimConfig::K_EOS = K_EOS_gui_;
        SimConfig::EPSILON_SQUARED = EPSILON_SQUARED_gui_;
        SimConfig::MAX_EXPECTED_PARTICLE_SPEED = MAX_EXPECTED_PARTICLE_SPEED_gui_;
        SimConfig::GAS_PARTICLE_JITTER_AMOUNT = GAS_PARTICLE_JITTER_AMOUNT_gui_;
        SimConfig::THETA_BARNES_HUT_SQUARED = static_cast<double>(theta_bh_gui_ * theta_bh_gui_);
    }
}



void Simulation::stepPhysics(double dt) {
    size_t numParticles = particle_list_.size();
    if (numParticles == 0) return;

    // Krok 0: Budowa struktur przyspieszających
    barnes_hut_tree_.build(particle_list_);
    sph_grid_object_.build(particle_list_, SimConfig::H_SPH);

    // Krok 1: Obliczanie gęstości i ciśnienia dla cząstek SPH
    if (sph_grid_object_.isBuilt()) { // Poprawione użycie isBuilt()
        for (size_t i = 0; i < numParticles; ++i) {
            if (particle_list_[i].type == ParticleType::GAS) {
                particle_list_[i].density = 0.0;
                const Particle& p_i = particle_list_[i];
                int cell_x_i, cell_y_i, cell_z_i;
                if (sph_grid_object_.getCellCoordinates(p_i.position, cell_x_i, cell_y_i, cell_z_i)) {
                    std::vector<size_t> neighbors = sph_grid_object_.getParticlesInNeighboringCells(cell_x_i, cell_y_i, cell_z_i);
                    for (size_t idx_j : neighbors) {
                        const Particle& p_j = particle_list_[idx_j];
                        double dist = (p_i.position - p_j.position).length();
                        particle_list_[i].density += p_j.mass * SPHUtils::kernelSPH(dist, SimConfig::H_SPH);
                    }
                }
                particle_list_[i].pressure = SPHUtils::calculatePressure(particle_list_[i].density, SimConfig::K_EOS, SimConfig::GAMMA_EOS);
            }
        }
    } else {
        for (size_t i = 0; i < numParticles; ++i) {
            if (particle_list_[i].type == ParticleType::GAS) {
                particle_list_[i].density = 0.0;
                particle_list_[i].pressure = 0.0;
            }
        }
    }
    
    // Krok 2: Obliczanie sił
    std::vector<Vector3d> forceList(numParticles, Vector3d(0.0, 0.0, 0.0));
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 1;
    if (numParticles < num_threads * 10) num_threads = 1;

    std::vector<std::thread> threads; // std::thread powinno być teraz rozpoznane
    threads.reserve(num_threads);
    size_t chunk_size = numParticles / num_threads;
    size_t remainder = numParticles % num_threads;
    size_t current_start_idx = 0;
    
    // Przekazujemy particle_list_ jako stałą referencję do wątków,
    // ponieważ odczytujemy z niej, a siły zapisujemy do osobnego forceList.
    const std::vector<Particle>& particleList_const_ref = particle_list_;

for (unsigned int i_thread = 0; i_thread < num_threads; ++i_thread) {
        size_t current_chunk_size = chunk_size + (i_thread < remainder ? 1 : 0);
        if (current_start_idx >= numParticles) break;
        size_t end_idx = std::min(current_start_idx + current_chunk_size, numParticles);
        
        threads.emplace_back( // Wywołanie metody klasy w nowym wątku
            &Simulation::calculateForcesForChunk, this, 
            current_start_idx,
            end_idx,
            std::ref(particleList_const_ref), 
            std::ref(forceList),            
            std::ref(barnes_hut_tree_),     
            std::ref(sph_grid_object_)      
        );
        current_start_idx = end_idx;
    }
    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }

    for (size_t i = 0; i < numParticles; ++i) {
        if (particle_list_[i].mass == 0.0) continue;
        Vector3d acceleration = forceList[i] / particle_list_[i].mass;
        particle_list_[i].velocity = particle_list_[i].velocity + (acceleration * dt);
        
        if (particle_list_[i].type == ParticleType::GAS) {
            particle_list_[i].trail_history.push_front(particle_list_[i].position);
            if (particle_list_[i].trail_history.size() > MAX_TRAIL_HISTORY_SIM) {
                particle_list_[i].trail_history.pop_back();
            }
        }
        particle_list_[i].position = particle_list_[i].position + (particle_list_[i].velocity * dt);
    }
}

// Metoda calculateForcesForChunk musi być metodą klasy Simulation lub przyjmować więcej parametrów
void Simulation::calculateForcesForChunk(
    size_t start_idx, size_t end_idx,
    const std::vector<Particle>& particles_read_only,
    std::vector<Vector3d>& forces_write,
    const Octree& tree, const SPHGrid& sph_grid_ref
) {
    for (size_t i = start_idx; i < end_idx; ++i) {
        forces_write[i] = Vector3d(0, 0, 0);
        const Particle& p_i = particles_read_only[i];

        // 1. Siły grawitacyjne (Barnes-Hut)
        forces_write[i] = forces_write[i] + tree.calculateForce(p_i, p_i.id);

        // 2. Siły ciśnienia SPH
        if (p_i.type == ParticleType::GAS && sph_grid_ref.isBuilt()) {
            int cell_x_i, cell_y_i, cell_z_i;
            if (sph_grid_ref.getCellCoordinates(p_i.position, cell_x_i, cell_y_i, cell_z_i)) {
                std::vector<size_t> neighbors = sph_grid_ref.getParticlesInNeighboringCells(cell_x_i, cell_y_i, cell_z_i);
                for (size_t idx_j : neighbors) {
                    if (p_i.id == particles_read_only[idx_j].id) continue;
                    const Particle& p_j = particles_read_only[idx_j];
                    forces_write[i] = forces_write[i] + SPHUtils::calculatePressureForceSymmetric(p_i, p_j, SimConfig::H_SPH, SimConfig::K_EOS, SimConfig::GAMMA_EOS);
                }
            }
        }
    }
}


void Simulation::renderFrame(float frame_delta_time_seconds) {
    if (!renderer_ptr_) return;
    renderer_ptr_->beginFrame();

    if (background_manager_ptr_) {
        background_manager_ptr_->draw(renderer_ptr_->getSDLRenderer(), camera_,
                                      renderer_ptr_->getWindowWidth(), 
                                      renderer_ptr_->getWindowHeight(),
                                      frame_delta_time_seconds);
    }
    renderer_ptr_->drawParticles(particle_list_, camera_, simulation_initial_size_);

    ImGui_ImplSDLRenderer2_NewFrame();
    ImGui_ImplSDL2_NewFrame();
    ImGui::NewFrame();
    setupGUI();
    ImGui::Render();
    ImGui_ImplSDLRenderer2_RenderDrawData(ImGui::GetDrawData(), renderer_ptr_->getSDLRenderer());

    renderer_ptr_->endFrame();
}

void Simulation::setupGUI() {
    ImGuiIO& io = ImGui::GetIO();
    ImGui::Begin("Simulation Control Panel");
    ImGui::Text("Particle Counts:");
    ImGui::InputInt("Dark Matter", &num_dark_matter_particles_gui_, 100, 1000);
    num_dark_matter_particles_gui_ = std::max(0, num_dark_matter_particles_gui_);
    ImGui::InputInt("Gas", &num_gas_particles_gui_, 100, 1000);
    num_gas_particles_gui_ = std::max(0, num_gas_particles_gui_);
    ImGui::Separator();
    ImGui::Text("Physics Parameters:");

    ImGui::SliderFloat("Time Step (dt)", &time_step_gui_, 0.0001f, 0.1f, "%.4f");
    ImGui::SliderFloat("Gravitational Constant (G)", &G_gui_, 0.1f, 10.0f); // Użyj pola klasy
    ImGui::SliderFloat("Grav. Softening (eps^2)", &EPSILON_SQUARED_gui_, 0.001f, 1.0f, "%.3f", ImGuiSliderFlags_Logarithmic); // Użyj pola klasy
    ImGui::SliderFloat("SPH Smoothing Length (H)", &H_SPH_gui_, 1.0f, 20.0f); // Użyj pola klasy
    ImGui::SliderFloat("EOS Constant (K)", &K_EOS_gui_, 0.01f, 1.0f); // Użyj pola klasy
    ImGui::SliderFloat("Barnes-Hut THETA", &theta_bh_gui_, 0.1f, 1.5f, "%.2f"); // Użyj pola klasy
    ImGui::SliderFloat("Max Expected Speed (Trails)", &MAX_EXPECTED_PARTICLE_SPEED_gui_, 0.1f, 20.0f); // Użyj pola klasy
    ImGui::SliderFloat("Gas Jitter Amount", &GAS_PARTICLE_JITTER_AMOUNT_gui_, 0.0f, 2.0f); // Użyj pola klasy

    ImGui::Separator();
    ImGui::Text("Camera:");
    ImGui::SliderFloat("Camera Smoothness", &camera_.lerp_speed, 0.01f, 0.3f);
    ImGui::Text("Position: (%.1f, %.1f, %.1f)", camera_.position.x, camera_.position.y, camera_.position.z);
    ImGui::Text("R: %.1f, Az: %.2f, El: %.2f", camera_.current_radius, camera_.current_azimuth, camera_.current_elevation);
    ImGui::Separator();
    if (ImGui::Button("Restart Simulation (R)")) {
        restartSimulation();
    }
    ImGui::SameLine();
    if (ImGui::Button(simulation_paused_ ? "Resume (P)" : "Pause (P)")) {
        simulation_paused_ = !simulation_paused_;
    }
    ImGui::Separator();
    ImGui::Text("Simulation Time: %.2f", current_time_);
    ImGui::Text("Frame: %u", main_loop_counter_); // Użyj pola klasy
    if (total_physics_steps_ > 0) {
        ImGui::Text("Avg. Physics Step Time: %.2f ms", static_cast<double>(total_physics_time_ms_) / total_physics_steps_); // Użyj pola klasy
    }
    ImGui::Text("FPS: %.1f", io.Framerate);
    ImGui::End();
}

void Simulation::initializeParticlesInternal(int num_dm, int num_gas, double initial_size, double total_mass) {
    particle_list_.clear();
    particle_list_.reserve(static_cast<size_t>(num_dm) + num_gas);
    if (num_dm + num_gas == 0) return;

    double totalParticlesDouble = static_cast<double>(num_dm + num_gas);
    double darkMatterMassFraction = (num_dm > 0) ? (static_cast<double>(num_dm) / totalParticlesDouble) : 0.0;
    double gasMassFraction = (num_gas > 0) ? (static_cast<double>(num_gas) / totalParticlesDouble) : 0.0;
    double massPerDarkMatterParticle = (num_dm > 0) ? (total_mass * darkMatterMassFraction) / num_dm : 0.0;
    double massPerGasParticle = (num_gas > 0) ? (total_mass * gasMassFraction) / num_gas : 0.0;
    
    size_t current_id = 0;
    for (int i = 0; i < num_dm; ++i) {
        particle_list_.emplace_back(current_id++, ParticleType::DARK_MATTER, massPerDarkMatterParticle,
            Vector3d(RandomUtils::generateRandomNumberFromTo(-initial_size / 2.0, initial_size / 2.0),
                     RandomUtils::generateRandomNumberFromTo(-initial_size / 2.0, initial_size / 2.0),
                     RandomUtils::generateRandomNumberFromTo(-initial_size / 2.0, initial_size / 2.0)),
            Vector3d(RandomUtils::generateRandomNumberFromTo(-0.000001 * initial_size, 0.000001 * initial_size),
                     RandomUtils::generateRandomNumberFromTo(-0.000001 * initial_size, 0.000001 * initial_size),
                     RandomUtils::generateRandomNumberFromTo(-0.000001 * initial_size, 0.000001 * initial_size))
        );
    }
    for (int i = 0; i < num_gas; ++i) {
         particle_list_.emplace_back(current_id++, ParticleType::GAS, massPerGasParticle,
            Vector3d(RandomUtils::generateRandomNumberFromTo(-initial_size / 2.0, initial_size / 2.0),
                     RandomUtils::generateRandomNumberFromTo(-initial_size / 2.0, initial_size / 2.0),
                     RandomUtils::generateRandomNumberFromTo(-initial_size / 2.0, initial_size / 2.0)),
            Vector3d(RandomUtils::generateRandomNumberFromTo(-0.000001 * initial_size, 0.000001 * initial_size),
                     RandomUtils::generateRandomNumberFromTo(-0.000001 * initial_size, 0.000001 * initial_size),
                     RandomUtils::generateRandomNumberFromTo(-0.000001 * initial_size, 0.000001 * initial_size))
        );
    }
}

void Simulation::restartSimulation() {
    initializeParticlesInternal(num_dark_matter_particles_gui_, num_gas_particles_gui_, simulation_initial_size_, 1e5);
    if (background_manager_ptr_) {
        background_manager_ptr_->initialize(simulation_initial_size_, 10.0, SimConfig::DUST_FIELD_SCALE);
    }
    current_time_ = 0.0;
    main_loop_counter_ = 0;
    total_physics_time_ms_ = 0;
    total_physics_steps_ = 0;
    camera_.initializeSpherical(initial_camera_radius_, initial_camera_azimuth_, initial_camera_elevation_);
    simulation_paused_ = false;
    // std::cout << "Info: Simulation restarted." << std::endl;
}