#pragma once

#include "Vector3d.h"

struct Camera {
    Vector3d position;
    Vector3d lookAt; 
    Vector3d up;     
    double fovY_degrees; 
    double nearClip;
    double farClip;

    double current_radius;
    double current_azimuth;
    double current_elevation;

    double target_radius;
    double target_azimuth;
    double target_elevation;

    float lerp_speed; // Camera smoothing speed


    Camera() : position(0, 0, 400), lookAt(0,0,0), up(0,1,0), fovY_degrees(60.0), nearClip(1.0), farClip(2000.0),
               current_radius(400.0), current_azimuth(0.0), current_elevation(0.0),
               target_radius(400.0), target_azimuth(0.0), target_elevation(0.0),
               lerp_speed(0.1f) 
               {}

    void initializeSpherical(double radius, double azimuth, double elevation);
    void updatePositionFromSpherical();
    void smoothUpdate();
};