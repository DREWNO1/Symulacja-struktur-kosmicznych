#include "Camera.h"
#include "SimConfig.h"

void Camera::initializeSpherical(double radius, double azimuth, double elevation)
{
    current_radius = target_radius = radius;
    current_azimuth = target_azimuth = azimuth;
    current_elevation = target_elevation = elevation;
    updatePositionFromSpherical();
}

void Camera::updatePositionFromSpherical()
{
    position.x = current_radius * std::cos(current_elevation) * std::sin(current_azimuth);
    position.y = current_radius * std::sin(current_elevation);
    position.z = current_radius * std::cos(current_elevation) * std::cos(current_azimuth);
}

void Camera::smoothUpdate()
{
    current_radius += (target_radius - current_radius) * lerp_speed;

    double delta_azimuth = target_azimuth - current_azimuth;
    while (delta_azimuth <= -SimConfig::PI)
        delta_azimuth += 2.0 * SimConfig::PI;
    while (delta_azimuth > SimConfig::PI)
        delta_azimuth -= 2.0 * SimConfig::PI;
    current_azimuth += delta_azimuth * lerp_speed;

    current_azimuth = std::fmod(current_azimuth, 2.0 * SimConfig::PI);
    if (current_azimuth < 0)
        current_azimuth += 2.0 * SimConfig::PI;

    current_elevation += (target_elevation - current_elevation) * lerp_speed;
    current_elevation = std::max(-SimConfig::PI / 2.0 + 0.01, std::min(SimConfig::PI / 2.0 - 0.01, current_elevation));

    updatePositionFromSpherical();
}