
#include "BackgroundManager.h" 
#include "Renderer.h"          
#include "RandomUtils.h" 
#include "SimConfig.h" 
#include <algorithm> 
#include <vector>    


Star::Star() 
    : position(0,0,0), base_brightness(0), size(1), 
      twinkle_timer(0.0f), current_brightness_offset(0.0f), twinkle_duration(0.0f) {}

void Star::initialize(double field_width, double field_height, double field_depth) {
    position.x = RandomUtils::generateRandomNumberFromTo(-field_width / 2.0, field_width / 2.0);
    position.y = RandomUtils::generateRandomNumberFromTo(-field_height / 2.0, field_height / 2.0);
    position.z = RandomUtils::generateRandomNumberFromTo(-field_depth / 2.0, field_depth / 2.0);
    base_brightness = static_cast<Uint8>(RandomUtils::generateRandomNumberFromTo(30, 150));
    double rand_size_roll = RandomUtils::generateRandomNumberFromTo(0.0, 1.0);
    if (rand_size_roll < 0.7) size = 1;
    else if (rand_size_roll < 0.9) size = 2;
    else size = 3;
    twinkle_timer = RandomUtils::generateRandomFloatFromTo(0.0f, 5.0f);
    current_brightness_offset = 0.0f;
    twinkle_duration = 0.0f;
}

void Star::updateTwinkle(float deltaTime) {
    twinkle_timer -= deltaTime;
    if (twinkle_timer <= 0.0f) {
        if (twinkle_duration <= 0.0f) { 
            current_brightness_offset = RandomUtils::generateRandomFloatFromTo(-60.0f, 60.0f);
            twinkle_duration = RandomUtils::generateRandomFloatFromTo(0.1f, 0.5f);
            twinkle_timer = twinkle_duration + RandomUtils::generateRandomFloatFromTo(1.0f, 10.0f); 
        } else { 
            current_brightness_offset = 0.0f;
            twinkle_duration = 0.0f; 
            twinkle_timer = RandomUtils::generateRandomFloatFromTo(1.0f, 10.0f); 
        }
    }

    if (twinkle_duration > 0.0f) { 
        twinkle_duration -= deltaTime;
        if (twinkle_duration <= 0.0f) {
            current_brightness_offset = 0.0f; 
        }
    }
}


CosmicDustParticle::CosmicDustParticle()
    : position(0,0,0), velocity(0,0,0), alpha(0), size(1) {}

void CosmicDustParticle::initialize(double dust_spread, double max_speed) {
    position.x = RandomUtils::generateRandomNumberFromTo(-dust_spread / 2.0, dust_spread / 2.0);
    position.y = RandomUtils::generateRandomNumberFromTo(-dust_spread / 2.0, dust_spread / 2.0);
    position.z = RandomUtils::generateRandomNumberFromTo(-dust_spread / 2.0, dust_spread / 2.0);
    velocity.x = RandomUtils::generateRandomNumberFromTo(-max_speed, max_speed);
    velocity.y = RandomUtils::generateRandomNumberFromTo(-max_speed, max_speed);
    velocity.z = RandomUtils::generateRandomNumberFromTo(-max_speed, max_speed);
    alpha = static_cast<Uint8>(RandomUtils::generateRandomNumberFromTo(10, 40));
    size = (RandomUtils::generateRandomNumberFromTo(0.0, 1.0) < 0.8) ? 1 : 2;
}

void CosmicDustParticle::updatePosition(double deltaTime, double bounds) {
    position = position + velocity * deltaTime; 
    if (position.x > bounds) position.x = -bounds + (position.x - bounds);
    if (position.x < -bounds) position.x = bounds - (-bounds - position.x);
    if (position.y > bounds) position.y = -bounds + (position.y - bounds);
    if (position.y < -bounds) position.y = bounds - (-bounds - position.y);
    if (position.z > bounds) position.z = -bounds + (position.z - bounds);
    if (position.z < -bounds) position.z = bounds - (-bounds - position.z);
}






BackgroundManager::BackgroundManager() {}

void BackgroundManager::initialize(double simulation_size, double starfield_depth_multiplier, double dust_field_scale_param) {
    starfield_.clear();
    starfield_.reserve(NUM_STARS);
    double star_field_width_height = simulation_size * starfield_depth_multiplier;
    double star_field_depth = SimConfig::STARFIELD_DEPTH; 

    for (int i = 0; i < NUM_STARS; ++i) {
        Star s;
        s.initialize(star_field_width_height, star_field_width_height, star_field_depth);
        starfield_.push_back(s);
    }

    cosmic_dust_layer_.clear();
    cosmic_dust_layer_.reserve(NUM_DUST_PARTICLES);
    double dust_spread = simulation_size * dust_field_scale_param; 
    double dust_max_speed = SimConfig::DUST_MAX_SPEED;

    for (int i = 0; i < NUM_DUST_PARTICLES; ++i) {
        CosmicDustParticle dp;
        dp.initialize(dust_spread, dust_max_speed);
        cosmic_dust_layer_.push_back(dp);
    }
}

void BackgroundManager::update(float frame_delta_time_seconds, double simulation_size, double dust_field_scale_param) {
    for (auto& star : starfield_) {
        star.updateTwinkle(frame_delta_time_seconds);
    }

    double dust_bounds = simulation_size * dust_field_scale_param / 2.0;
    for (auto& dp : cosmic_dust_layer_) {
        dp.updatePosition(frame_delta_time_seconds, dust_bounds);
    }
}


void BackgroundManager::draw(Renderer& renderer_obj, const Camera& camera, float frame_delta_time_seconds) const {
    SDL_Renderer* sdl_renderer = renderer_obj.getSDLRenderer(); 
    if (!sdl_renderer) return;

    int screen_width = renderer_obj.getWindowWidth();
    int screen_height = renderer_obj.getWindowHeight();

    SDL_SetRenderDrawBlendMode(sdl_renderer, SDL_BLENDMODE_BLEND);

    
    for (const auto& star : starfield_) {
        
        ProjectedPoint pp_star = renderer_obj.project3DTo2D(star.position, camera); 
        if (pp_star.visible) {
            int final_brightness_val = static_cast<int>(star.base_brightness + star.current_brightness_offset);
            Uint8 final_brightness = static_cast<Uint8>(std::max(0, std::min(255, final_brightness_val)));
            SDL_SetRenderDrawColor(sdl_renderer, final_brightness, final_brightness, final_brightness, 255);
            if (star.size == 1) {
                SDL_RenderDrawPoint(sdl_renderer, pp_star.screenX, pp_star.screenY);
            } else {
                SDL_Rect star_rect = {pp_star.screenX - star.size / 2, pp_star.screenY - star.size / 2, star.size, star.size};
                SDL_RenderFillRect(sdl_renderer, &star_rect);
            }
        }
    }

    
    for (const auto& dp : cosmic_dust_layer_) {
        ProjectedPoint pp_dust = renderer_obj.project3DTo2D(dp.position, camera);
        if (pp_dust.visible) {
            SDL_SetRenderDrawColor(sdl_renderer, 50, 50, 60, dp.alpha);
            if (dp.size == 1) {
                SDL_RenderDrawPoint(sdl_renderer, pp_dust.screenX, pp_dust.screenY);
            } else {
                SDL_Rect dust_rect = {pp_dust.screenX - dp.size / 2, pp_dust.screenY - dp.size / 2, dp.size, dp.size};
                SDL_RenderFillRect(sdl_renderer, &dust_rect);
            }
        }
    }
}