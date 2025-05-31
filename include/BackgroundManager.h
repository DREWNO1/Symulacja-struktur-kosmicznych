// include/BackgroundManager.h
#ifndef BACKGROUNDMANAGER_H
#define BACKGROUNDMANAGER_H

#include <vector>
#include <SDL2/SDL_stdinc.h> // Dla Uint8
#include "Vector3d.h"
#include "Camera.h"

class Renderer; // Forward declaration - to jest poprawne dla pliku .h

struct Star {
    Vector3d position;
    Uint8 base_brightness;
    int size;
    float twinkle_timer;
    float current_brightness_offset;
    float twinkle_duration;

    Star();
    void initialize(double field_width, double field_height, double field_depth);
    void updateTwinkle(float deltaTime);
};

struct CosmicDustParticle {
    Vector3d position;
    Vector3d velocity;
    Uint8 alpha;
    int size;

    CosmicDustParticle();
    void initialize(double dust_spread, double max_speed);
    void updatePosition(double deltaTime, double bounds);
};

class BackgroundManager {
public:
    BackgroundManager();

    void initialize(double simulation_size, double starfield_depth_multiplier, double dust_field_scale);
    void update(float frame_delta_time_seconds, double simulation_size, double dust_field_scale);
    
    // Deklaracja metody draw
    void draw(Renderer& renderer_obj, const Camera& camera, float frame_delta_time_seconds) const;

private:
    std::vector<Star> starfield_;
    std::vector<CosmicDustParticle> cosmic_dust_layer_;

    static const int NUM_STARS = 500; 
    static const int NUM_DUST_PARTICLES = 1000;
};

#endif // BACKGROUNDMANAGER_H