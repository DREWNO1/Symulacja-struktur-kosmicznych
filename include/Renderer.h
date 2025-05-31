
#ifndef RENDERER_H
#define RENDERER_H

#include <vector>
#include <string>
#include <SDL2/SDL.h>
#include "Vector3d.h" 
#include "Camera.h"
#include "Particle.h"
#include "SimConfig.h"

struct ProjectedPoint {
    int screenX;
    int screenY;
    double depth;
    bool visible;
};

class Renderer {
public:
    Renderer();
    ~Renderer();

    bool initialize(const char* title, int window_width, int window_height);
    void shutdown();

    void beginFrame();
    void drawParticles(const std::vector<Particle>& particleList, const Camera& camera, double simulation_initial_size);
    void endFrame();

    SDL_Renderer* getSDLRenderer() const { return sdl_renderer_; }
    SDL_Window* getSDLWindow() const { return sdl_window_; }
    int getWindowWidth() const { return window_width_; }
    int getWindowHeight() const { return window_height_; }

    ProjectedPoint project3DTo2D(const Vector3d& world_pos, const Camera& camera) const;

private:
    SDL_Window* sdl_window_;
    SDL_Renderer* sdl_renderer_;
    int window_width_;
    int window_height_;

    
    
};

#endif 