#include <iostream>     // For std::cout, std::cerr, std::endl
#include <vector>       // For std::vector
#include <random>       // For std::mt19937, std::uniform_real_distribution
#include <cmath>        // For std::sqrt, std::pow, std::sin, std::cos, std::tan, std::max, std::min
#include <limits>       // For std::numeric_limits
#include <utility>      // For std::pair
#include <algorithm>    // For std::min, std::max
#include <string>       // For std::string
#include <thread>       // For std::thread, std::thread::hardware_concurrency
#include <deque>        // For std::deque (particle trail history)
#include <memory>       // For std::unique_ptr

// --- SDL2 Headers ---
#include <SDL2/SDL.h>   // Main SDL library

// --- ImGui Headers ---
#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_sdlrenderer2.h"

// --- Custom Simulation Headers ---
#include "SimConfig.h"   // Simulation parameters
#include "Vector3d.h"    // 3D vector operations
#include "Camera.h"      // Camera and projection
#include "Particle.h"    // Particle structure definition
#include "Octree.h"      // Octree class definition
#include "SPHGrid.h"

// --- Global Random Number Engine Helper ---
// (Przeniesiono z poprzedniego main.cpp, jeśli nie ma go w innym wspólnym pliku utility)
std::mt19937& get_random_engine() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

// --- Camera and Projection Structures (ProjectedPoint pozostaje, bo jest specyficzne dla logiki renderowania w main) ---
struct ProjectedPoint {
    int screenX;
    int screenY;
    double depth;
    bool visible;
};

// --- Particle System Structures (ParticleType i Particle są teraz w Particle.h) ---
// Enum ParticleType i struct Particle są teraz zdefiniowane w Particle.h
// MAX_TRAIL_HISTORY może pozostać jako globalna stała tutaj lub zostać przeniesiona do SimConfig.h
// albo stać się częścią konfiguracji cząstki lub funkcji rysującej.
// W Particle.h jest jako komentarz, więc dla spójności z poprzednim kodem, zostawmy go tutaj.
const int MAX_TRAIL_HISTORY = 10; //

// --- Star structure for background (with twinkling) ---
struct Star {
    Vector3d position;
    Uint8 base_brightness;
    int size;
    float twinkle_timer;
    float current_brightness_offset;
    float twinkle_duration;
};

std::vector<Star> starfield; //

const int NUM_STARS = 500; //
const double STARFIELD_DEPTH = 1500.0; //

// --- Cosmic Dust structure (with velocity) ---
struct CosmicDustParticle {
    Vector3d position;
    Vector3d velocity;
    Uint8 alpha;
    int size;
};

std::vector<CosmicDustParticle> cosmicDustLayer; //

const int NUM_DUST_PARTICLES = 1000; //
const double DUST_FIELD_SCALE = 1.5; //
const double DUST_MAX_SPEED = 0.05; //

// --- Utility Functions (te, które nie są częścią klas) ---
double generateRandomNumberFromTo(double min_val, double max_val) { //
    if (min_val > max_val) std::swap(min_val, max_val); //
    std::uniform_real_distribution<double> distrib(min_val, max_val); //
    return distrib(get_random_engine()); //
}

float generateRandomFloatFromTo(float min_val, float max_val) { //
    if (min_val > max_val) std::swap(min_val, max_val); //
    std::uniform_real_distribution<float> distrib(min_val, max_val); //
    return distrib(get_random_engine()); //
}

// calculateDistanceSquared, calculateDistance, calculateUnitVector mogą być używane przez metody Vector3d
// (np. (v1-v2).lengthSquared()), więc mogą nie być już potrzebne jako wolne funkcje,
// jeśli Vector3d.h/cpp zapewnia pełną funkcjonalność.
// Oryginalny main.cpp miał je, ale Vector3d.cpp już implementuje lengthSquared(), length(), normalized().
// calculateGravityForce została przeniesiona do octree.cpp jako calculateDirectGravityForce.

// --- SPH Kernel Functions (pozostają globalne lub mogą być w namespace SPHUtils) ---
double kernelSPH(double r, double h) { //
    if (r < 0.0 || r > h) { //
        return 0.0; //
    }
    if (h == 0.0) return 0.0; //

    double h_sq = h * h; //
    double r_sq = r * r; //
    double h_ninth = h_sq * h_sq * h_sq * h_sq * h; //
    if (h_ninth == 0.0) return 0.0; //

    double factor = 315.0 / (64.0 * SimConfig::PI * h_ninth); //
    double term = h_sq - r_sq; //
    return factor * term * term * term; //
}

double derivativeKernelSPH(double r, double h) { //
    if (r < 0.0 || r > h) { //
        return 0.0; //
    }
    if (h == 0.0) return 0.0; //

    double h_sixth = std::pow(h, 6.0); //
    if (h_sixth == 0.0) return 0.0; //

    double factor = -45.0 / (SimConfig::PI * h_sixth); //
    double term = h - r; //
    return factor * term * term; //
}

// --- Barnes-Hut Octree ---
// Definicje OctreeNode i Octree są teraz w Octree.h i Octree.cpp
// Globalna instancja drzewa
Octree barnes_hut_tree; //

SPHGrid sph_grid_object;

// --- Initialization Functions ---
void initializeStarfield(double field_width, double field_height, double field_depth) { //
    starfield.clear(); //
    starfield.reserve(NUM_STARS); //
    for (int i = 0; i < NUM_STARS; ++i) { //
        Star s; //
        s.position.x = generateRandomNumberFromTo(-field_width / 2.0, field_width / 2.0); //
        s.position.y = generateRandomNumberFromTo(-field_height / 2.0, field_height / 2.0); //
        s.position.z = generateRandomNumberFromTo(-field_depth / 2.0, field_depth / 2.0); //
        s.base_brightness = static_cast<Uint8>(generateRandomNumberFromTo(30, 150)); //
        double rand_size_roll = generateRandomNumberFromTo(0.0, 1.0); //
        if (rand_size_roll < 0.7) s.size = 1; //
        else if (rand_size_roll < 0.9) s.size = 2; //
        else s.size = 3; //
        s.twinkle_timer = generateRandomFloatFromTo(0.0f, 5.0f); //
        s.current_brightness_offset = 0.0f; //
        s.twinkle_duration = 0.0f; //
        starfield.push_back(s); //
    }
}

void initializeCosmicDust(double simulation_size) { //
    cosmicDustLayer.clear(); //
    cosmicDustLayer.reserve(NUM_DUST_PARTICLES); //
    double dust_spread = simulation_size * DUST_FIELD_SCALE; //
    for (int i = 0; i < NUM_DUST_PARTICLES; ++i) { //
        CosmicDustParticle dp; //
        dp.position.x = generateRandomNumberFromTo(-dust_spread / 2.0, dust_spread / 2.0); //
        dp.position.y = generateRandomNumberFromTo(-dust_spread / 2.0, dust_spread / 2.0); //
        dp.position.z = generateRandomNumberFromTo(-dust_spread / 2.0, dust_spread / 2.0); //
        dp.velocity.x = generateRandomNumberFromTo(-DUST_MAX_SPEED, DUST_MAX_SPEED); //
        dp.velocity.y = generateRandomNumberFromTo(-DUST_MAX_SPEED, DUST_MAX_SPEED); //
        dp.velocity.z = generateRandomNumberFromTo(-DUST_MAX_SPEED, DUST_MAX_SPEED); //
        dp.alpha = static_cast<Uint8>(generateRandomNumberFromTo(10, 40)); //
        dp.size = (generateRandomNumberFromTo(0.0, 1.0) < 0.8) ? 1 : 2; //
        cosmicDustLayer.push_back(dp); //
    }
}

void updateCosmicDust(double deltaTime, double simulation_size) { //
    double bounds = simulation_size * DUST_FIELD_SCALE / 2.0; //
    for (auto& dp : cosmicDustLayer) { //
        dp.position = dp.position + dp.velocity * deltaTime; //
        if (dp.position.x > bounds) dp.position.x = -bounds + (dp.position.x - bounds); //
        if (dp.position.x < -bounds) dp.position.x = bounds - (-bounds - dp.position.x); //
        if (dp.position.y > bounds) dp.position.y = -bounds + (dp.position.y - bounds); //
        if (dp.position.y < -bounds) dp.position.y = bounds - (-bounds - dp.position.y); //
        if (dp.position.z > bounds) dp.position.z = -bounds + (dp.position.z - bounds); //
        if (dp.position.z < -bounds) dp.position.z = bounds - (-bounds - dp.position.z); //
    }
}

std::vector<Particle> initializeParticles(int p_numDarkMatter, int p_numGas, double p_initialSize, double p_totalMass) { //
    std::vector<Particle> particleList; //
    particleList.reserve(static_cast<size_t>(p_numDarkMatter) + p_numGas); //
    if(p_numDarkMatter + p_numGas == 0) return particleList; //

    double totalParticlesDouble = static_cast<double>(p_numDarkMatter + p_numGas); //
    double darkMatterMassFraction = (p_numDarkMatter > 0) ? (static_cast<double>(p_numDarkMatter) / totalParticlesDouble) : 0.0; //
    double gasMassFraction = (p_numGas > 0) ? (static_cast<double>(p_numGas) / totalParticlesDouble) : 0.0; //
    double massPerDarkMatterParticle = (p_numDarkMatter > 0) ? (p_totalMass * darkMatterMassFraction) / p_numDarkMatter : 0.0; //
    double massPerGasParticle = (p_numGas > 0) ? (p_totalMass * gasMassFraction) / p_numGas : 0.0; //
    
    size_t current_id = 0; //

    // Inicjalizacja cząstek ciemnej materii
    for (int i = 0; i < p_numDarkMatter; ++i) { //
        // Używamy konstruktora z parametrami z Particle.cpp/Particle.h
        Particle newParticle(current_id++, ParticleType::DARK_MATTER, massPerDarkMatterParticle); //
        newParticle.position.x = generateRandomNumberFromTo(-p_initialSize/2.0, p_initialSize/2.0); //
        newParticle.position.y = generateRandomNumberFromTo(-p_initialSize/2.0, p_initialSize/2.0); //
        newParticle.position.z = generateRandomNumberFromTo(-p_initialSize/2.0, p_initialSize/2.0); //
        newParticle.velocity.x = generateRandomNumberFromTo(-0.000001 * p_initialSize, 0.000001 * p_initialSize); //
        newParticle.velocity.y = generateRandomNumberFromTo(-0.000001 * p_initialSize, 0.000001 * p_initialSize); //
        newParticle.velocity.z = generateRandomNumberFromTo(-0.000001 * p_initialSize, 0.000001 * p_initialSize); //
        // density i pressure są już 0.0 z konstruktora
        particleList.push_back(newParticle); //
    }
    // Inicjalizacja cząstek gazu
    for (int i = 0; i < p_numGas; ++i) { //
        Particle newParticle(current_id++, ParticleType::GAS, massPerGasParticle); //
        newParticle.position.x = generateRandomNumberFromTo(-p_initialSize/2.0, p_initialSize/2.0); //
        newParticle.position.y = generateRandomNumberFromTo(-p_initialSize/2.0, p_initialSize/2.0); //
        newParticle.position.z = generateRandomNumberFromTo(-p_initialSize/2.0, p_initialSize/2.0); //
        newParticle.velocity.x = generateRandomNumberFromTo(-0.000001 * p_initialSize, 0.000001 * p_initialSize); //
        newParticle.velocity.y = generateRandomNumberFromTo(-0.000001 * p_initialSize, 0.000001 * p_initialSize); //
        newParticle.velocity.z = generateRandomNumberFromTo(-0.000001 * p_initialSize, 0.000001 * p_initialSize); //
        particleList.push_back(newParticle); //
    }
    return particleList; //
}


// --- Physics Calculation Functions ---
double calculatePressureSPH(double density) { //
    if (density < 0.0) density = 0.0; //
    return SimConfig::K_EOS * std::pow(density, SimConfig::GAMMA_EOS); //
}

Vector3d calculatePressureForceSPHSymmetric(const Particle& p_p1, const Particle& p_p2) { //
    // double distance = calculateDistance(p_p1.position, p_p2.position); // Użycie metody Vector3d
    double distance = (p_p1.position - p_p2.position).length(); //
    if (distance >= 2.0 * SimConfig::H_SPH || distance == 0.0) return Vector3d(0.0, 0.0, 0.0); //

    Vector3d positionVector1ToVector2 = p_p2.position - p_p1.position; //
    // Vector3d unitVector1ToVector2 = calculateUnitVector(positionVector1ToVector2); // Użycie metody Vector3d
    Vector3d unitVector1ToVector2 = positionVector1ToVector2.normalized(); //
    double pressure1 = p_p1.pressure; //
    double pressure2 = p_p2.pressure; //
    double rho1 = p_p1.density; //
    double rho2 = p_p2.density; //
    if (rho1 == 0.0 || rho2 == 0.0) return Vector3d(0.0, 0.0, 0.0); //
    double pressureTerm = (pressure1 / (rho1 * rho1)) + (pressure2 / (rho2 * rho2)); //
    double deriv_W = derivativeKernelSPH(distance, SimConfig::H_SPH); //
    double scalarFactor = -p_p2.mass * pressureTerm * deriv_W; //
    return unitVector1ToVector2 * scalarFactor; //
}

void calculateForcesForChunk(
    size_t start_idx,
    size_t end_idx,
    const std::vector<Particle>& particles_read_only,
    std::vector<Vector3d>& forces_write,
    const Octree& tree,
    const SPHGrid& sph_grid_ref // Nowy argument
) {
    for (size_t i = start_idx; i < end_idx; ++i) { //
        forces_write[i] = Vector3d(0,0,0); //
        const Particle& p_i = particles_read_only[i]; //

        // 1. Siły grawitacyjne (Barnes-Hut)
        forces_write[i] = forces_write[i] + tree.calculateForce(p_i, p_i.id); //


        // 2. Siły ciśnienia SPH (z wykorzystaniem siatki SPH)
        if (p_i.type == ParticleType::GAS && sph_grid_ref.isBuilt()) { // Użyj obiektu siatki
            int cell_x_i, cell_y_i, cell_z_i;
            if (sph_grid_ref.getCellCoordinates(p_i.position, cell_x_i, cell_y_i, cell_z_i)) {
                std::vector<size_t> potential_neighbors_indices =
                    sph_grid_ref.getParticlesInNeighboringCells(cell_x_i, cell_y_i, cell_z_i);

                for (size_t particle_idx_j : potential_neighbors_indices) {
                    if (p_i.id == particles_read_only[particle_idx_j].id) continue;
                    const Particle& p_j = particles_read_only[particle_idx_j];
                    Vector3d pressureForce = calculatePressureForceSPHSymmetric(p_i, p_j);
                    forces_write[i] = forces_write[i] + pressureForce;
                }
            }
        }
    }
}


// --- Main Simulation Step ---
void simulateStep(std::vector<Particle>& particleList, double deltaTime) { //
    size_t numParticles = particleList.size(); //
    if (numParticles == 0) return; //

    barnes_hut_tree.build(particleList); //
    sph_grid_object.build(particleList, SimConfig::H_SPH);

if (sph_grid_object.isBuilt()) { // Sprawdź, czy siatka została poprawnie zbudowana
    for (size_t i = 0; i < numParticles; ++i) {
        if (particleList[i].type == ParticleType::GAS) {
            particleList[i].density = 0.0;
            const Particle& p_i = particleList[i];
            int cell_x_i, cell_y_i, cell_z_i;

            // Użyj metody klasy SPHGrid do uzyskania współrzędnych komórki
            if (sph_grid_object.getCellCoordinates(p_i.position, cell_x_i, cell_y_i, cell_z_i)) {
                // Pobierz cząstki z sąsiednich komórek (własnej i 26 otaczających)
                std::vector<size_t> potential_neighbors_indices = 
                    sph_grid_object.getParticlesInNeighboringCells(cell_x_i, cell_y_i, cell_z_i);

                for (size_t particle_idx_j : potential_neighbors_indices) {
                    const Particle& p_j = particleList[particle_idx_j];
                    double distance = (p_i.position - p_j.position).length();
                    particleList[i].density += p_j.mass * kernelSPH(distance, SimConfig::H_SPH);
                }
            }
            particleList[i].pressure = calculatePressureSPH(particleList[i].density);
        }
    }
    } else {
        for (size_t i = 0; i < numParticles; ++i) { // Ten blok 'else' może pozostać, jeśli siatka nie jest zbudowana
            if (particleList[i].type == ParticleType::GAS) {
                particleList[i].density = 0.0;
                particleList[i].pressure = 0.0;
            }
        }
    }

    std::vector<Vector3d> forceList(numParticles, Vector3d(0.0, 0.0, 0.0)); //
    unsigned int num_threads = std::thread::hardware_concurrency(); //
    if (num_threads == 0) num_threads = 1; //
    if (numParticles < num_threads * 10) num_threads = 1; //

    std::vector<std::thread> threads; //
    threads.reserve(num_threads); //
    size_t chunk_size = numParticles / num_threads; //
    size_t remainder = numParticles % num_threads; //
    size_t current_start_idx = 0; //
    const std::vector<Particle>& particleList_const_ref = particleList; //

    for (unsigned int i_thread = 0; i_thread < num_threads; ++i_thread) { //
        size_t current_chunk_size = chunk_size + (i_thread < remainder ? 1 : 0); //
        if (current_start_idx >= numParticles) break; //
        size_t end_idx = std::min(current_start_idx + current_chunk_size, numParticles); //
        threads.emplace_back(
            calculateForcesForChunk,
            current_start_idx,
            end_idx,
            std::ref(particleList_const_ref),
            std::ref(forceList),
            std::ref(barnes_hut_tree),
            std::ref(sph_grid_object) // Przekaż referencję do obiektu siatki
        );
        current_start_idx = end_idx; //
    }
    for (auto& t : threads) { //
        if (t.joinable()) { //
            t.join(); //
        }
    }

    for (size_t i = 0; i < numParticles; ++i) { //
        if (particleList[i].mass == 0.0) continue; //
        Vector3d acceleration = forceList[i] / particleList[i].mass; //
        particleList[i].velocity = particleList[i].velocity + (acceleration * deltaTime); //
        if (particleList[i].type == ParticleType::GAS) { //
            particleList[i].trail_history.push_front(particleList[i].position); //
            if (particleList[i].trail_history.size() > MAX_TRAIL_HISTORY) { //
                particleList[i].trail_history.pop_back(); //
            }
        }
        particleList[i].position = particleList[i].position + (particleList[i].velocity * deltaTime); //
    }
}

// --- SDL and Rendering Functions ---
std::pair<SDL_Window*, SDL_Renderer*> initializeSDL(int p_w, int p_h, const char* title) { //
    if (SDL_Init(SDL_INIT_VIDEO) < 0) { //
        std::cerr << "ERROR: SDL_Init failed: " << SDL_GetError() << std::endl; //
        return {nullptr, nullptr}; //
    }
    SDL_Window* window = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, p_w, p_h, SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE); //
    if (window == nullptr) { //
        std::cerr << "ERROR: SDL_CreateWindow failed: " << SDL_GetError() << std::endl; //
        SDL_Quit(); //
        return {nullptr, nullptr}; //
    }
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC); //
    if (renderer == nullptr) { //
        std::cerr << "ERROR: SDL_CreateRenderer failed: " << SDL_GetError() << std::endl; //
        SDL_DestroyWindow(window); //
        SDL_Quit(); //
        return {nullptr, nullptr}; //
    }
    return {window, renderer}; //
}

void closeSDL(SDL_Window* window, SDL_Renderer* renderer) { //
    if (renderer != nullptr) SDL_DestroyRenderer(renderer); //
    if (window != nullptr) SDL_DestroyWindow(window); //
    SDL_Quit(); //
}

Vector3d crossProduct(const Vector3d& a, const Vector3d& b) { //
    return Vector3d(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); //
}
double dotProduct(const Vector3d& a, const Vector3d& b) { //
    return a.x * b.x + a.y * b.y + a.z * b.z; //
}

ProjectedPoint project3DTo2D(const Vector3d& world_pos, const Camera& camera, int screen_width, int screen_height) { //
    ProjectedPoint p{}; //
    p.visible = false; //
    Vector3d D = (camera.lookAt - camera.position).normalized(); //
    Vector3d R = crossProduct(camera.up, D).normalized(); //
    Vector3d U = crossProduct(D, R); //
    Vector3d rel_pos = world_pos - camera.position; //
    Vector3d view_coords(dotProduct(rel_pos, R), dotProduct(rel_pos, U), dotProduct(rel_pos, D)); //

    if (view_coords.z < camera.nearClip || view_coords.z > camera.farClip) return p; //
    if (view_coords.z == 0) return p; //

    double fov_rad_y = camera.fovY_degrees * SimConfig::PI / 180.0; //
    double aspectRatio = static_cast<double>(screen_width) / screen_height; //
    double sy_proj = 1.0 / std::tan(fov_rad_y / 2.0); //
    double sx_proj = sy_proj / aspectRatio; //
    double proj_x_ndc = (view_coords.x * sx_proj) / view_coords.z; //
    double proj_y_ndc = (view_coords.y * sy_proj) / view_coords.z; //

    p.screenX = static_cast<int>((proj_x_ndc + 1.0) * 0.5 * screen_width); //
    p.screenY = static_cast<int>((1.0 - proj_y_ndc) * 0.5 * screen_height); //
    p.depth = view_coords.z; //
    p.visible = true; //
    return p; //
}

void drawStarfieldSDL(SDL_Renderer* renderer, const Camera& camera, int screen_width, int screen_height, float deltaTime) { //
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND); //
    for (auto& star : starfield) { //
        star.twinkle_timer -= deltaTime; //
        if (star.twinkle_timer <= 0.0f) { //
            if (star.twinkle_duration <= 0.0f) { //
                star.current_brightness_offset = generateRandomFloatFromTo(-60.0f, 60.0f); //
                star.twinkle_duration = generateRandomFloatFromTo(0.1f, 0.5f); //
                star.twinkle_timer = star.twinkle_duration + generateRandomFloatFromTo(1.0f, 10.0f); //
            } else {
                star.current_brightness_offset = 0.0f; //
                star.twinkle_duration = 0.0f; //
                star.twinkle_timer = generateRandomFloatFromTo(1.0f, 10.0f); //
            }
        }
        if (star.twinkle_duration > 0.0f) { //
            star.twinkle_duration -= deltaTime; //
            if (star.twinkle_duration <= 0.0f) { //
                 star.current_brightness_offset = 0.0f; //
            }
        }
        ProjectedPoint pp_star = project3DTo2D(star.position, camera, screen_width, screen_height); //
        if (pp_star.visible) { //
            int final_brightness_val = static_cast<int>(star.base_brightness + star.current_brightness_offset); //
            Uint8 final_brightness = static_cast<Uint8>(std::max(0, std::min(255, final_brightness_val))); //
            SDL_SetRenderDrawColor(renderer, final_brightness, final_brightness, final_brightness, 255); //
            if (star.size == 1) { //
                SDL_RenderDrawPoint(renderer, pp_star.screenX, pp_star.screenY); //
            } else {
                SDL_Rect star_rect = {pp_star.screenX - star.size / 2, pp_star.screenY - star.size / 2, star.size, star.size}; //
                SDL_RenderFillRect(renderer, &star_rect); //
            }
        }
    }
}

void drawCosmicDustSDL(SDL_Renderer* renderer, const Camera& camera, int screen_width, int screen_height) { //
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND); //
    for (const auto& dp : cosmicDustLayer) { //
        ProjectedPoint pp_dust = project3DTo2D(dp.position, camera, screen_width, screen_height); //
        if (pp_dust.visible) { //
            SDL_SetRenderDrawColor(renderer, 50, 50, 60, dp.alpha); //
            if (dp.size == 1) { //
                SDL_RenderDrawPoint(renderer, pp_dust.screenX, pp_dust.screenY); //
            } else {
                SDL_Rect dust_rect = {pp_dust.screenX - dp.size / 2, pp_dust.screenY - dp.size / 2, dp.size, dp.size}; //
                SDL_RenderFillRect(renderer, &dust_rect); //
            }
        }
    }
}

void drawParticlesSDL(SDL_Renderer* renderer, const std::vector<Particle>& particleList, //
                      const Camera& camera, int screen_width, int screen_height, double simulationInitialSize) { //
    double pressure_color_factor = 0.5; //
    SDL_BlendMode original_blend_mode; //
    SDL_GetRenderDrawBlendMode(renderer, &original_blend_mode); //

    for (const auto& particle : particleList) { //
        ProjectedPoint pp = project3DTo2D(particle.position, camera, screen_width, screen_height); //
        int jitter_x = 0, jitter_y = 0; //
        if (particle.type == ParticleType::GAS) { //
            jitter_x = static_cast<int>(generateRandomFloatFromTo(-SimConfig::GAS_PARTICLE_JITTER_AMOUNT, SimConfig::GAS_PARTICLE_JITTER_AMOUNT)); //
            jitter_y = static_cast<int>(generateRandomFloatFromTo(-SimConfig::GAS_PARTICLE_JITTER_AMOUNT, SimConfig::GAS_PARTICLE_JITTER_AMOUNT)); //
        }
        pp.screenX += jitter_x; //
        pp.screenY += jitter_y; //
        if (!pp.visible) continue; //

        double base_particle_size_dm = 3.0; //
        double base_particle_size_gas = 25.0; //
        double reference_depth = simulationInitialSize; //
        int particle_render_size; //
        if (particle.type == ParticleType::GAS) { //
            particle_render_size = static_cast<int>(std::max(2.0, base_particle_size_gas * reference_depth / (pp.depth + reference_depth*0.5) )); //
            particle_render_size = std::min(particle_render_size, 35); //
        } else {
            particle_render_size = static_cast<int>(std::max(1.0, base_particle_size_dm * reference_depth / (pp.depth + reference_depth*0.5) )); //
            particle_render_size = std::min(particle_render_size, 5); //
        }
        particle_render_size = std::max(1, particle_render_size); //
        Uint8 r_col_particle = 0, g_col_particle = 0, b_col_particle = 0; //

        if (particle.type == ParticleType::DARK_MATTER) { //
            r_col_particle = 100; g_col_particle = 100; b_col_particle = 120; //
            SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND); //
            SDL_SetRenderDrawColor(renderer, r_col_particle, g_col_particle, b_col_particle, 180); //
            if (particle_render_size <= 2) { //
                SDL_RenderDrawPoint(renderer, pp.screenX, pp.screenY); //
                if (particle_render_size == 2) { //
                    SDL_RenderDrawPoint(renderer, pp.screenX+1, pp.screenY); //
                    SDL_RenderDrawPoint(renderer, pp.screenX, pp.screenY+1); //
                    SDL_RenderDrawPoint(renderer, pp.screenX+1, pp.screenY+1); //
                }
            } else {
                SDL_Rect rect = {pp.screenX - particle_render_size / 2, pp.screenY - particle_render_size / 2, particle_render_size, particle_render_size}; //
                SDL_RenderFillRect(renderer, &rect); //
            }
        } else if (particle.type == ParticleType::GAS) { //
            double normalized_pressure = particle.pressure * pressure_color_factor; //
            if (normalized_pressure < 0.0) normalized_pressure = 0.0; //
            if (normalized_pressure < 0.2) { //
                double t = normalized_pressure / 0.2; //
                r_col_particle = static_cast<Uint8>(50 + t * 50); //
                g_col_particle = static_cast<Uint8>(50 + t * 50); //
                b_col_particle = static_cast<Uint8>(150 + t * 105); //
            } else if (normalized_pressure < 0.5) { //
                double t = (normalized_pressure - 0.2) / 0.3; //
                r_col_particle = static_cast<Uint8>(100 + t * 155); //
                g_col_particle = static_cast<Uint8>(100 - t * 50); //
                b_col_particle = static_cast<Uint8>(255 - t * 100); //
            } else if (normalized_pressure < 0.8) { //
                double t = (normalized_pressure - 0.5) / 0.3; //
                r_col_particle = 255; //
                g_col_particle = static_cast<Uint8>(50 + t * 205); //
                b_col_particle = static_cast<Uint8>(155 - t * 155); //
            } else {
                double t = std::min(1.0, (normalized_pressure - 0.8) / 0.2); //
                r_col_particle = 255; //
                g_col_particle = 255; //
                b_col_particle = static_cast<Uint8>(std::min(255.0, 100.0 + t * 155.0)); //
            }
            r_col_particle = std::min(static_cast<Uint8>(255), std::max(static_cast<Uint8>(0), r_col_particle)); //
            g_col_particle = std::min(static_cast<Uint8>(255), std::max(static_cast<Uint8>(0), g_col_particle)); //
            b_col_particle = std::min(static_cast<Uint8>(255), std::max(static_cast<Uint8>(0), b_col_particle)); //

            SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND); //
            int trail_idx = 0; //
            // double speed = particle.velocity.length(); // Użycie metody Vector3d
            double speed = particle.velocity.length(); //
            float speed_factor = std::min(1.0f, static_cast<float>(speed / SimConfig::MAX_EXPECTED_PARTICLE_SPEED)); //
            int num_trail_segments_to_draw = static_cast<int>(particle.trail_history.size() * (0.5f + speed_factor * 0.5f) ); //
            num_trail_segments_to_draw = std::min(num_trail_segments_to_draw, (int)particle.trail_history.size()); //

            for (const auto& trail_pos : particle.trail_history) { //
                if (trail_idx >= num_trail_segments_to_draw) break; //
                ProjectedPoint pp_trail = project3DTo2D(trail_pos, camera, screen_width, screen_height); //
                if (pp_trail.visible) { //
                    float age_factor = static_cast<float>(trail_idx) / MAX_TRAIL_HISTORY; //
                    Uint8 trail_alpha = static_cast<Uint8>(100 * (1.0f - age_factor) * (0.6f + speed_factor * 0.4f) ); //
                    int trail_render_size = std::max(1, static_cast<int>(particle_render_size * 0.3f * (1.0f - age_factor * 0.7f) )); //
                    Uint8 tr_r = static_cast<Uint8>(r_col_particle * (1.0f - age_factor * 0.6f) + 30 * age_factor * 0.6f); //
                    Uint8 tr_g = static_cast<Uint8>(g_col_particle * (1.0f - age_factor * 0.8f) + 30 * age_factor * 0.8f); //
                    Uint8 tr_b = static_cast<Uint8>(b_col_particle * (1.0f - age_factor * 0.4f) + 80 * age_factor * 0.4f); //
                    SDL_SetRenderDrawColor(renderer, tr_r, tr_g, tr_b, trail_alpha); //
                    SDL_Rect trail_rect = {pp_trail.screenX - trail_render_size / 2, pp_trail.screenY - trail_render_size / 2, trail_render_size, trail_render_size}; //
                    SDL_RenderFillRect(renderer, &trail_rect); //
                }
                trail_idx++; //
            }

            SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND); //
            float density_alpha_factor = std::min(1.0f, 0.1f + static_cast<float>(particle.density / (SimConfig::H_SPH * 5.0) ) ); //
            Uint8 gas_alpha = static_cast<Uint8>(40 * density_alpha_factor); //
            gas_alpha = std::max((Uint8)10, std::min((Uint8)100, gas_alpha)); //
            int num_gas_layers = 2; //
            for(int i_gas_layer = 0; i_gas_layer < num_gas_layers; ++i_gas_layer) { //
                float size_factor = 1.0f - (float)i_gas_layer / (num_gas_layers + 1); //
                int current_gas_size = std::max(1, (int)(particle_render_size * size_factor * 1.2f)); //
                Uint8 current_gas_alpha = static_cast<Uint8>(gas_alpha * (1.0f - (float)i_gas_layer / (num_gas_layers * 2.0f))); //
                current_gas_alpha = std::max((Uint8)5, current_gas_alpha); //
                SDL_SetRenderDrawColor(renderer, r_col_particle, g_col_particle, b_col_particle, current_gas_alpha); //
                SDL_Rect gas_rect = {pp.screenX - current_gas_size / 2, pp.screenY - current_gas_size / 2, current_gas_size, current_gas_size}; //
                SDL_RenderFillRect(renderer, &gas_rect); //
            }

            SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_ADD); //
            double glow_intensity_factor = 1.0 + std::min(1.5, std::max(0.0, normalized_pressure - 0.6) * 2.5); //
            int num_glow_layers = 2; //
            double glow_size_multiplier_base = 0.7 * glow_intensity_factor; //
            double glow_alpha_base = 35.0 * glow_intensity_factor; //
            glow_alpha_base = std::min(100.0, glow_alpha_base); //
            for (int i_glow = 0; i_glow < num_glow_layers; ++i_glow) { //
                double current_glow_multiplier = glow_size_multiplier_base + i_glow * 0.4 * glow_intensity_factor; //
                int glow_size = static_cast<int>(particle_render_size * current_glow_multiplier); //
                glow_size = std::max(particle_render_size / 2, glow_size); //
                Uint8 glow_alpha = static_cast<Uint8>(glow_alpha_base / (i_glow + 1)); //
                glow_alpha = std::max((Uint8)10, glow_alpha); //
                SDL_SetRenderDrawColor(renderer, r_col_particle, g_col_particle, b_col_particle, glow_alpha); //
                SDL_Rect glow_rect = {pp.screenX - glow_size / 2, pp.screenY - glow_size / 2, glow_size, glow_size}; //
                SDL_RenderFillRect(renderer, &glow_rect); //
            }

            if (normalized_pressure > 0.9) { //
                SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_ADD); //
                int core_render_size = std::max(1, particle_render_size / 3); //
                core_render_size = std::min(core_render_size, 4); //
                SDL_SetRenderDrawColor(renderer, 255, 255, 255, 220); //
                SDL_Rect core_rect = {pp.screenX - core_render_size / 2, pp.screenY - core_render_size / 2, core_render_size, core_render_size}; //
                SDL_RenderFillRect(renderer, &core_rect); //
            }
        }
    }
    SDL_SetRenderDrawBlendMode(renderer, original_blend_mode); //
}


// --- Main Function ---
#ifdef main //
#undef main //
#endif
int main(int argc, char* argv[]) { //
    int window_width = 1920; //
    int window_height = 1080; //
    auto sdl_components = initializeSDL(window_width, window_height, "Cosmic Structure Simulation 3D - Volumetric Gas"); //
    SDL_Window* window = sdl_components.first; //
    SDL_Renderer* renderer = sdl_components.second; //
    if (window == nullptr || renderer == nullptr) { //
        std::cerr << "FATAL ERROR: SDL Initialization failed. Exiting." << std::endl; //
        return 1; //
    }

    IMGUI_CHECKVERSION(); //
    ImGui::CreateContext(); //
    ImGuiIO& io = ImGui::GetIO(); (void)io; //
    ImGui::StyleColorsDark(); //
    ImGui_ImplSDL2_InitForSDLRenderer(window, renderer); //
    ImGui_ImplSDLRenderer2_Init(renderer); //

    Camera camera; //
    double simulationInitialSize = 200.0; //
    double initial_camera_radius = simulationInitialSize * 2.25; //
    double initial_camera_azimuth = 0.0; //
    double initial_camera_elevation = 0.0; //
    camera.initializeSpherical(initial_camera_radius, initial_camera_azimuth, initial_camera_elevation); //
    camera.lookAt = Vector3d(0, 0, 0); //
    camera.up = Vector3d(0, 1, 0); //
    camera.fovY_degrees = 60.0; //
    camera.nearClip = simulationInitialSize / 100.0; //
    camera.farClip = simulationInitialSize * 10.0; //
    camera.lerp_speed = 0.08f; //

    initializeStarfield(simulationInitialSize * 10, simulationInitialSize * 10, STARFIELD_DEPTH * 2.0); //
    initializeCosmicDust(simulationInitialSize); //

    bool mouse_left_button_down = false; //
    int mouse_prev_x = 0, mouse_prev_y = 0; //
    const double MOUSE_SENSITIVITY = 0.005; //
    const double MOUSE_WHEEL_SENSITIVITY = simulationInitialSize / 20.0; //

    static int numDarkMatterParticles_gui = 500; //
    static int numGasParticles_gui = 2000; //
    static float timeStep_gui = 0.005f; //
    static float G_gui = static_cast<float>(SimConfig::G); //
    static float H_SPH_gui = static_cast<float>(SimConfig::H_SPH); //
    static float K_EOS_gui = static_cast<float>(SimConfig::K_EOS); //
    static float EPSILON_SQUARED_gui = static_cast<float>(SimConfig::EPSILON_SQUARED); //
    static float MAX_EXPECTED_PARTICLE_SPEED_gui = SimConfig::MAX_EXPECTED_PARTICLE_SPEED; //
    static float GAS_PARTICLE_JITTER_AMOUNT_gui = SimConfig::GAS_PARTICLE_JITTER_AMOUNT; //
    // THETA_BARNES_HUT_SQUARED jest teraz w SimConfig.h, więc inicjalizacja THETA_BH_gui jest poprawna.
    static float THETA_BH_gui = static_cast<float>(std::sqrt(SimConfig::THETA_BARNES_HUT_SQUARED)); //

    std::vector<Particle> particleList = initializeParticles(numDarkMatterParticles_gui, numGasParticles_gui, simulationInitialSize, 1e5); //

    double currentTime = 0.0; //
    bool running = true; //
    bool simulation_paused = false; //
    Uint32 mainLoopCounter = 0; //
    Uint32 total_physics_time_ms = 0; //
    int total_physics_steps = 0; //
    Uint32 last_frame_ticks = SDL_GetTicks(); //

    std::cout << "Info: Simulation started (3D with ImGui, Volumetric Gas, Dynamic SPH Grid, Barnes-Hut Gravity)." << std::endl; //
    std::cout << "Info: Controls: Arrow keys (orbit), A/Z (zoom), Mouse LMB+drag (orbit), Mouse Wheel (zoom), R (reset), P (pause)" << std::endl; //

    while (running) { //
        Uint32 current_ticks = SDL_GetTicks(); //
        float frame_delta_time_seconds = (current_ticks - last_frame_ticks) / 1000.0f; //
        if (frame_delta_time_seconds <= 0) frame_delta_time_seconds = 1.0f/60.0f; //
        if (frame_delta_time_seconds > 0.1f) frame_delta_time_seconds = 0.1f; //
        last_frame_ticks = current_ticks; //

        SDL_Event event; //
        while (SDL_PollEvent(&event)) { //
            ImGui_ImplSDL2_ProcessEvent(&event); //
            bool imgui_captures_mouse = io.WantCaptureMouse; //
            bool imgui_captures_keyboard = io.WantCaptureKeyboard; //
            if (event.type == SDL_QUIT) running = false; //
            if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_RESIZED) { //
                SDL_GetWindowSize(window, &window_width, &window_height); //
            }
            if (!imgui_captures_mouse) { //
                if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_LEFT) { //
                    mouse_left_button_down = true; //
                    SDL_GetMouseState(&mouse_prev_x, &mouse_prev_y); //
                }
                if (event.type == SDL_MOUSEBUTTONUP && event.button.button == SDL_BUTTON_LEFT) { //
                    mouse_left_button_down = false; //
                }
                if (event.type == SDL_MOUSEMOTION && mouse_left_button_down) { //
                    int current_mouse_x, current_mouse_y; //
                    SDL_GetMouseState(&current_mouse_x, &current_mouse_y); //
                    int delta_x = current_mouse_x - mouse_prev_x; //
                    int delta_y = current_mouse_y - mouse_prev_y; //
                    camera.target_azimuth -= static_cast<double>(delta_x) * MOUSE_SENSITIVITY; //
                    camera.target_elevation -= static_cast<double>(delta_y) * MOUSE_SENSITIVITY; //
                    mouse_prev_x = current_mouse_x; //
                    mouse_prev_y = current_mouse_y; //
                }
                if (event.type == SDL_MOUSEWHEEL) { //
                    if (event.wheel.y > 0) camera.target_radius -= MOUSE_WHEEL_SENSITIVITY; //
                    else if (event.wheel.y < 0) camera.target_radius += MOUSE_WHEEL_SENSITIVITY; //
                }
            }
            if (!imgui_captures_keyboard) { //
                if (event.type == SDL_KEYDOWN) { //
                    if (event.key.keysym.sym == SDLK_r) { //
                        particleList = initializeParticles(numDarkMatterParticles_gui, numGasParticles_gui, simulationInitialSize, 1e5); //
                        initializeStarfield(simulationInitialSize * 10, simulationInitialSize * 10, STARFIELD_DEPTH * 2.0); //
                        initializeCosmicDust(simulationInitialSize); //
                        currentTime = 0.0; mainLoopCounter = 0; total_physics_time_ms = 0; total_physics_steps = 0; //
                        camera.initializeSpherical(initial_camera_radius, initial_camera_azimuth, initial_camera_elevation); //
                        std::cout << "Info: Simulation restarted." << std::endl; //
                    }
                    if (event.key.keysym.sym == SDLK_p) { //
                        simulation_paused = !simulation_paused; //
                        std::cout << "Info: Simulation " << (simulation_paused ? "paused." : "resumed.") << std::endl; //
                    }
                    if (event.key.keysym.sym == SDLK_LEFT) camera.target_azimuth -= camera.lerp_speed * 5.0f; //
                    if (event.key.keysym.sym == SDLK_RIGHT) camera.target_azimuth += camera.lerp_speed * 5.0f; //
                    if (event.key.keysym.sym == SDLK_UP) camera.target_elevation += camera.lerp_speed * 5.0f; //
                    if (event.key.keysym.sym == SDLK_DOWN) camera.target_elevation -= camera.lerp_speed * 5.0f; //
                    if (event.key.keysym.sym == SDLK_a) camera.target_radius -= MOUSE_WHEEL_SENSITIVITY; //
                    if (event.key.keysym.sym == SDLK_z) camera.target_radius += MOUSE_WHEEL_SENSITIVITY; //
                }
            }
            camera.target_elevation = std::max(-SimConfig::PI / 2.0 + 0.02, std::min(SimConfig::PI / 2.0 - 0.02, camera.target_elevation)); //
            camera.target_radius = std::max(camera.nearClip * 3.0, camera.target_radius); //
        }

        camera.smoothUpdate(); //
        if (!simulation_paused) { //
            updateCosmicDust(frame_delta_time_seconds, simulationInitialSize); //
        }

        if (timeStep_gui > 0 && !simulation_paused) { //
            SimConfig::G = G_gui; //
            SimConfig::H_SPH = H_SPH_gui; //
            SimConfig::K_EOS = K_EOS_gui; //
            SimConfig::EPSILON_SQUARED = EPSILON_SQUARED_gui; //
            SimConfig::MAX_EXPECTED_PARTICLE_SPEED = MAX_EXPECTED_PARTICLE_SPEED_gui; //
            SimConfig::GAS_PARTICLE_JITTER_AMOUNT = GAS_PARTICLE_JITTER_AMOUNT_gui; //
            SimConfig::THETA_BARNES_HUT_SQUARED = static_cast<double>(THETA_BH_gui * THETA_BH_gui); //


            Uint32 sim_step_start_ticks = SDL_GetTicks(); //
            simulateStep(particleList, static_cast<double>(timeStep_gui)); //
            Uint32 sim_step_duration_ms = SDL_GetTicks() - sim_step_start_ticks; //
            total_physics_time_ms += sim_step_duration_ms; //
            total_physics_steps++; //
            currentTime += static_cast<double>(timeStep_gui); //
        }
        mainLoopCounter++; //

        ImGui_ImplSDLRenderer2_NewFrame(); //
        ImGui_ImplSDL2_NewFrame(); //
        ImGui::NewFrame(); //
        ImGui::Begin("Simulation Control Panel"); //
        ImGui::Text("Particle Counts:"); //
        ImGui::InputInt("Dark Matter", &numDarkMatterParticles_gui, 100, 1000); //
        numDarkMatterParticles_gui = std::max(0, numDarkMatterParticles_gui); //
        ImGui::InputInt("Gas", &numGasParticles_gui, 100, 1000); //
        numGasParticles_gui = std::max(0, numGasParticles_gui); //
        ImGui::Separator(); //
        ImGui::Text("Physics Parameters:"); //
        ImGui::SliderFloat("Time Step (dt)", &timeStep_gui, 0.0001f, 0.1f, "%.4f"); //
        ImGui::SliderFloat("Gravitational Constant (G)", &G_gui, 0.1f, 10.0f); //
        ImGui::SliderFloat("Grav. Softening (eps^2)", &EPSILON_SQUARED_gui, 0.001f, 1.0f, "%.3f", ImGuiSliderFlags_Logarithmic); //
        if (ImGui::SliderFloat("SPH Smoothing Length (H)", &H_SPH_gui, 1.0f, 20.0f)) { //
             SimConfig::H_SPH = H_SPH_gui; //
        }
        ImGui::SliderFloat("EOS Constant (K)", &K_EOS_gui, 0.01f, 1.0f); //
        ImGui::SliderFloat("Barnes-Hut THETA", &THETA_BH_gui, 0.1f, 1.5f, "%.2f"); //
        ImGui::SliderFloat("Max Expected Speed (Trails)", &MAX_EXPECTED_PARTICLE_SPEED_gui, 0.1f, 20.0f); //
        ImGui::SliderFloat("Gas Jitter Amount", &GAS_PARTICLE_JITTER_AMOUNT_gui, 0.0f, 2.0f); //
        ImGui::Separator(); //
        ImGui::Text("Camera:"); //
        ImGui::SliderFloat("Camera Smoothness", &camera.lerp_speed, 0.01f, 0.3f); //
        ImGui::Text("Position: (%.1f, %.1f, %.1f)", camera.position.x, camera.position.y, camera.position.z); //
        ImGui::Text("R: %.1f, Az: %.2f, El: %.2f", camera.current_radius, camera.current_azimuth, camera.current_elevation); //
        ImGui::Separator(); //
        if (ImGui::Button("Restart Simulation (R)")) { //
            particleList = initializeParticles(numDarkMatterParticles_gui, numGasParticles_gui, simulationInitialSize, 1e5); //
            initializeStarfield(simulationInitialSize * 10, simulationInitialSize * 10, STARFIELD_DEPTH * 2.0); //
            initializeCosmicDust(simulationInitialSize); //
            currentTime = 0.0; mainLoopCounter = 0; total_physics_time_ms = 0; total_physics_steps = 0; //
            camera.initializeSpherical(initial_camera_radius, initial_camera_azimuth, initial_camera_elevation); //
            simulation_paused = false; //
            std::cout << "Info: Simulation restarted via ImGui." << std::endl; //
        }
        ImGui::SameLine(); //
        if (ImGui::Button(simulation_paused ? "Resume (P)" : "Pause (P)")) { //
            simulation_paused = !simulation_paused; //
        }
        ImGui::Separator(); //
        ImGui::Text("Simulation Time: %.2f", currentTime); //
        ImGui::Text("Frame: %u", mainLoopCounter); //
        if (total_physics_steps > 0) { //
            ImGui::Text("Avg. Physics Step Time: %.2f ms", static_cast<double>(total_physics_time_ms) / total_physics_steps); //
        }
        ImGui::Text("FPS: %.1f", io.Framerate); //
        ImGui::End(); //

        SDL_SetRenderDrawColor(renderer, 0, 0, 5, 255); //
        SDL_RenderClear(renderer); //
        drawStarfieldSDL(renderer, camera, window_width, window_height, frame_delta_time_seconds); //
        drawCosmicDustSDL(renderer, camera, window_width, window_height); //
        SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND); //
        drawParticlesSDL(renderer, particleList, camera, window_width, window_height, simulationInitialSize); //
        ImGui::Render(); //
        ImGui_ImplSDLRenderer2_RenderDrawData(ImGui::GetDrawData(), renderer); //
        SDL_RenderPresent(renderer); //

        Uint32 frame_render_duration_ms = SDL_GetTicks() - current_ticks; //
        const Uint32 TARGET_FRAME_DURATION_MS = 1000 / 60; //
        if (frame_render_duration_ms < TARGET_FRAME_DURATION_MS && !io.WantCaptureMouse && !io.WantCaptureKeyboard) { //
            SDL_Delay(TARGET_FRAME_DURATION_MS - frame_render_duration_ms); //
        }
    }

    ImGui_ImplSDLRenderer2_Shutdown(); //
    ImGui_ImplSDL2_Shutdown(); //
    ImGui::DestroyContext(); //
    closeSDL(window, renderer); //
    std::cout << "Simulation finished." << std::endl; //
    return 0; //
}