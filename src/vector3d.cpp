#include "Vector3d.h"
#include <stdexcept>

Vector3d Vector3d::operator/(double scalar) const {
    if (scalar == 0) {
        throw std::domain_error("Division by zero in Vector3d operator/");
    }
    return Vector3d(x / scalar, y / scalar, z / scalar);
}

double Vector3d::lengthSquared() const {
    return x * x + y * y + z * z;
}

double Vector3d::length() const {
    return std::sqrt(lengthSquared());
}

Vector3d Vector3d::normalized() const {
    double l = length();
    if (l == 0) {
        return Vector3d(0, 0, 0);
    }
    return Vector3d(x / l, y / l, z / l);
}