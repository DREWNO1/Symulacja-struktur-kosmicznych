
#pragma once

#include <iostream>
#include <cmath> 

struct Vector3d {
    double x;
    double y;
    double z;

    Vector3d(double p_x = 0.0, double p_y = 0.0, double p_z = 0.0)
        : x(p_x), y(p_y), z(p_z) {}

    
    Vector3d operator-(const Vector3d& other) const {
        return Vector3d(x - other.x, y - other.y, z - other.z);
    }

    Vector3d operator+(const Vector3d& other) const {
        return Vector3d(x + other.x, y + other.y, z + other.z);
    }

    Vector3d operator*(double scalar) const {
        return Vector3d(x * scalar, y * scalar, z * scalar);
    }


    Vector3d operator/(double scalar) const; 

    double lengthSquared() const; 
    double length() const;        
    Vector3d normalized() const;  

    
    double dot(const Vector3d& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vector3d cross(const Vector3d& other) const {
        return Vector3d(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }
};