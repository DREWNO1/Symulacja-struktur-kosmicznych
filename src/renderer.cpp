// src/renderer.cpp
#include "Renderer.h"
#include "SimConfig.h" // Dla SimConfig::PI
#include "RandomUtils.h" 
#include <algorithm>    
#include <cmath>        
#include <SDL2/SDL_log.h>

Renderer::Renderer()
    : sdl_window_(nullptr), sdl_renderer_(nullptr), window_width_(0), window_height_(0) {}

Renderer::~Renderer() {
    shutdown();
}

bool Renderer::initialize(const char* title, int window_width, int window_height) {
    window_width_ = window_width;
    window_height_ = window_height;

    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_Init failed: %s", SDL_GetError());
        return false;
    }

    sdl_window_ = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                 window_width_, window_height_, SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    if (sdl_window_ == nullptr) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_CreateWindow failed: %s", SDL_GetError());
        SDL_Quit();
        return false;
    }

    sdl_renderer_ = SDL_CreateRenderer(sdl_window_, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (sdl_renderer_ == nullptr) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "SDL_CreateRenderer failed: %s", SDL_GetError());
        SDL_DestroyWindow(sdl_window_);
        SDL_Quit();
        return false;
    }
    return true;
}

void Renderer::shutdown() {
    if (sdl_renderer_ != nullptr) {
        SDL_DestroyRenderer(sdl_renderer_);
        sdl_renderer_ = nullptr;
    }
    if (sdl_window_ != nullptr) {
        SDL_DestroyWindow(sdl_window_);
        sdl_window_ = nullptr;
    }
    // SDL_Quit(); // Powinno być wywołane raz, na końcu aplikacji
}

void Renderer::beginFrame() {
    if (!sdl_renderer_) return;
    SDL_GetWindowSize(sdl_window_, &window_width_, &window_height_);

    SDL_SetRenderDrawColor(sdl_renderer_, 0, 0, 5, 255);
    SDL_RenderClear(sdl_renderer_);
}

void Renderer::endFrame() {
    if (!sdl_renderer_) return;
    SDL_RenderPresent(sdl_renderer_);
}

// Usunięto definicje Renderer::crossProduct i Renderer::dotProduct

ProjectedPoint Renderer::project3DTo2D(const Vector3d& world_pos, const Camera& camera) const {
    ProjectedPoint p{};
    p.visible = false;

    Vector3d D = (camera.lookAt - camera.position).normalized();
    // Użyj metod z Vector3d
    Vector3d R = camera.up.cross(D).normalized(); 
    Vector3d U = D.cross(R);                     

    Vector3d rel_pos = world_pos - camera.position;
    // Użyj metod z Vector3d
    Vector3d view_coords(rel_pos.dot(R), rel_pos.dot(U), rel_pos.dot(D)); 

    if (view_coords.z < camera.nearClip || view_coords.z > camera.farClip) return p;
    if (view_coords.z == 0) return p;

    double fov_rad_y = camera.fovY_degrees * SimConfig::PI / 180.0;
    double aspect_ratio = (window_height_ > 0) ? (static_cast<double>(window_width_) / window_height_) : 1.0;
    double sy_proj = 1.0 / std::tan(fov_rad_y / 2.0);
    double sx_proj = sy_proj / aspect_ratio;

    double proj_x_ndc = (view_coords.x * sx_proj) / view_coords.z;
    double proj_y_ndc = (view_coords.y * sy_proj) / view_coords.z;

    p.screenX = static_cast<int>((proj_x_ndc + 1.0) * 0.5 * window_width_);
    p.screenY = static_cast<int>((1.0 - proj_y_ndc) * 0.5 * window_height_);
    p.depth = view_coords.z;
    p.visible = true;
    return p;
}

// Implementacja drawParticles (pełna, jak w poprzedniej odpowiedzi)
void Renderer::drawParticles(const std::vector<Particle>& particleList, const Camera& camera, double simulation_initial_size) {
    if (!sdl_renderer_) return;

    double pressure_color_factor = 0.5; 
    SDL_BlendMode original_blend_mode;
    SDL_GetRenderDrawBlendMode(sdl_renderer_, &original_blend_mode);

    for (const auto& particle : particleList) {
        ProjectedPoint pp = this->project3DTo2D(particle.position, camera); 
        
        int jitter_x = 0, jitter_y = 0;
        if (particle.type == ParticleType::GAS) {
            jitter_x = static_cast<int>(RandomUtils::generateRandomFloatFromTo(-SimConfig::GAS_PARTICLE_JITTER_AMOUNT, SimConfig::GAS_PARTICLE_JITTER_AMOUNT));
            jitter_y = static_cast<int>(RandomUtils::generateRandomFloatFromTo(-SimConfig::GAS_PARTICLE_JITTER_AMOUNT, SimConfig::GAS_PARTICLE_JITTER_AMOUNT));
        }
        pp.screenX += jitter_x;
        pp.screenY += jitter_y;

        if (!pp.visible) continue;

        double base_particle_size_dm = 3.0;
        double base_particle_size_gas = 25.0;
        double reference_depth = simulation_initial_size;
        int particle_render_size;

        if (particle.type == ParticleType::GAS) {
            particle_render_size = static_cast<int>(std::max(2.0, base_particle_size_gas * reference_depth / (pp.depth + reference_depth * 0.5)));
            particle_render_size = std::min(particle_render_size, 35);
        } else {
            particle_render_size = static_cast<int>(std::max(1.0, base_particle_size_dm * reference_depth / (pp.depth + reference_depth * 0.5)));
            particle_render_size = std::min(particle_render_size, 5);
        }
        particle_render_size = std::max(1, particle_render_size);
        Uint8 r_col_particle = 0, g_col_particle = 0, b_col_particle = 0;

        if (particle.type == ParticleType::DARK_MATTER) {
            r_col_particle = 100; g_col_particle = 100; b_col_particle = 120;
            SDL_SetRenderDrawBlendMode(sdl_renderer_, SDL_BLENDMODE_BLEND);
            SDL_SetRenderDrawColor(sdl_renderer_, r_col_particle, g_col_particle, b_col_particle, 180);
            if (particle_render_size <= 2) {
                SDL_RenderDrawPoint(sdl_renderer_, pp.screenX, pp.screenY);
                if (particle_render_size == 2) { 
                    SDL_RenderDrawPoint(sdl_renderer_, pp.screenX + 1, pp.screenY);
                    SDL_RenderDrawPoint(sdl_renderer_, pp.screenX, pp.screenY + 1);
                    SDL_RenderDrawPoint(sdl_renderer_, pp.screenX + 1, pp.screenY + 1);
                }
            } else {
                SDL_Rect rect = {pp.screenX - particle_render_size / 2, pp.screenY - particle_render_size / 2, particle_render_size, particle_render_size};
                SDL_RenderFillRect(sdl_renderer_, &rect);
            }
        } else if (particle.type == ParticleType::GAS) {
            double normalized_pressure = particle.pressure * pressure_color_factor;
            if (normalized_pressure < 0.0) normalized_pressure = 0.0;
            if (normalized_pressure < 0.2) {
                double t = normalized_pressure / 0.2;
                r_col_particle = static_cast<Uint8>(50 + t * 50); g_col_particle = static_cast<Uint8>(50 + t * 50); b_col_particle = static_cast<Uint8>(150 + t * 105);
            } else if (normalized_pressure < 0.5) {
                double t = (normalized_pressure - 0.2) / 0.3;
                r_col_particle = static_cast<Uint8>(100 + t * 155); g_col_particle = static_cast<Uint8>(100 - t * 50); b_col_particle = static_cast<Uint8>(255 - t * 100);
            } else if (normalized_pressure < 0.8) {
                double t = (normalized_pressure - 0.5) / 0.3;
                r_col_particle = 255; g_col_particle = static_cast<Uint8>(50 + t * 205); b_col_particle = static_cast<Uint8>(155 - t * 155);
            } else {
                double t = std::min(1.0, (normalized_pressure - 0.8) / 0.2);
                r_col_particle = 255; g_col_particle = 255; b_col_particle = static_cast<Uint8>(std::min(255.0, 100.0 + t * 155.0));
            }
            
            SDL_SetRenderDrawBlendMode(sdl_renderer_, SDL_BLENDMODE_BLEND);
            int trail_idx = 0;
            double speed = particle.velocity.length();
            float speed_factor = std::min(1.0f, static_cast<float>(speed / SimConfig::MAX_EXPECTED_PARTICLE_SPEED));
            int num_trail_segments_to_draw = static_cast<int>(particle.trail_history.size() * (0.5f + speed_factor * 0.5f));
            num_trail_segments_to_draw = std::min(num_trail_segments_to_draw, (int)particle.trail_history.size());
            
            const int MAX_TRAIL_HISTORY_CONST = 10; // Zastąp stałą z SimConfig jeśli tam jest, lub użyj lokalnej

            for (const auto& trail_pos : particle.trail_history) {
                if (trail_idx >= num_trail_segments_to_draw) break;
                ProjectedPoint pp_trail = this->project3DTo2D(trail_pos, camera);
                if (pp_trail.visible) {
                    float age_factor = static_cast<float>(trail_idx) / MAX_TRAIL_HISTORY_CONST;
                    Uint8 trail_alpha = static_cast<Uint8>(100 * (1.0f - age_factor) * (0.6f + speed_factor * 0.4f) );
                    int trail_render_size = std::max(1, static_cast<int>(particle_render_size * 0.3f * (1.0f - age_factor * 0.7f) ));
                    Uint8 tr_r = static_cast<Uint8>(r_col_particle * (1.0f - age_factor * 0.6f) + 30 * age_factor * 0.6f);
                    Uint8 tr_g = static_cast<Uint8>(g_col_particle * (1.0f - age_factor * 0.8f) + 30 * age_factor * 0.8f);
                    Uint8 tr_b = static_cast<Uint8>(b_col_particle * (1.0f - age_factor * 0.4f) + 80 * age_factor * 0.4f);
                    SDL_SetRenderDrawColor(sdl_renderer_, tr_r, tr_g, tr_b, trail_alpha);
                    SDL_Rect trail_rect = {pp_trail.screenX - trail_render_size / 2, pp_trail.screenY - trail_render_size / 2, trail_render_size, trail_render_size};
                    SDL_RenderFillRect(sdl_renderer_, &trail_rect);
                }
                trail_idx++;
            }

            SDL_SetRenderDrawBlendMode(sdl_renderer_, SDL_BLENDMODE_BLEND);
            float density_alpha_factor = std::min(1.0f, 0.1f + static_cast<float>(particle.density / (SimConfig::H_SPH * 5.0)));
            Uint8 gas_alpha = static_cast<Uint8>(40 * density_alpha_factor);
            gas_alpha = std::max((Uint8)10, std::min((Uint8)100, gas_alpha));
            int num_gas_layers = 2;
            for (int i_gas_layer = 0; i_gas_layer < num_gas_layers; ++i_gas_layer) {
                float size_factor = 1.0f - (float)i_gas_layer / (num_gas_layers + 1);
                int current_gas_size = std::max(1, (int)(particle_render_size * size_factor * 1.2f));
                Uint8 current_gas_alpha = static_cast<Uint8>(gas_alpha * (1.0f - (float)i_gas_layer / (num_gas_layers * 2.0f)));
                 current_gas_alpha = std::max((Uint8)5, current_gas_alpha);
                SDL_SetRenderDrawColor(sdl_renderer_, r_col_particle, g_col_particle, b_col_particle, current_gas_alpha);
                SDL_Rect gas_rect = {pp.screenX - current_gas_size / 2, pp.screenY - current_gas_size / 2, current_gas_size, current_gas_size};
                SDL_RenderFillRect(sdl_renderer_, &gas_rect);
            }

            SDL_SetRenderDrawBlendMode(sdl_renderer_, SDL_BLENDMODE_ADD);
            double glow_intensity_factor = 1.0 + std::min(1.5, std::max(0.0, normalized_pressure - 0.6) * 2.5);
            int num_glow_layers = 2;
            double glow_size_multiplier_base = 0.7 * glow_intensity_factor;
            double glow_alpha_base = 35.0 * glow_intensity_factor;
            glow_alpha_base = std::min(100.0, glow_alpha_base);
            for (int i_glow = 0; i_glow < num_glow_layers; ++i_glow) {
                double current_glow_multiplier = glow_size_multiplier_base + i_glow * 0.4 * glow_intensity_factor;
                int glow_size = static_cast<int>(particle_render_size * current_glow_multiplier);
                glow_size = std::max(particle_render_size / 2, glow_size);
                Uint8 glow_alpha = static_cast<Uint8>(glow_alpha_base / (i_glow + 1));
                glow_alpha = std::max((Uint8)10, glow_alpha);
                SDL_SetRenderDrawColor(sdl_renderer_, r_col_particle, g_col_particle, b_col_particle, glow_alpha);
                SDL_Rect glow_rect = {pp.screenX - glow_size / 2, pp.screenY - glow_size / 2, glow_size, glow_size};
                SDL_RenderFillRect(sdl_renderer_, &glow_rect);
            }

            if (normalized_pressure > 0.9) {
                SDL_SetRenderDrawBlendMode(sdl_renderer_, SDL_BLENDMODE_ADD);
                int core_render_size = std::max(1, particle_render_size / 3);
                core_render_size = std::min(core_render_size, 4);
                SDL_SetRenderDrawColor(sdl_renderer_, 255, 255, 255, 220);
                SDL_Rect core_rect = {pp.screenX - core_render_size / 2, pp.screenY - core_render_size / 2, core_render_size, core_render_size};
                SDL_RenderFillRect(sdl_renderer_, &core_rect);
            }
        }
    }
    SDL_SetRenderDrawBlendMode(sdl_renderer_, original_blend_mode);
}