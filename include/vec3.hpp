#ifndef VEC3_HPP
#define VEC3_HPP

#include <cmath>

/**
 * Lightweight 3D vector struct for HPC simulations.
 * Uses inline functions for performance-critical operations.
 */
struct Vec3 {
    double x, y, z;

    // Default constructor
    constexpr Vec3() : x(0.0), y(0.0), z(0.0) {}

    // Parameterized constructor
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // Addition
    inline Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    inline Vec3& operator+=(const Vec3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    // Subtraction
    inline Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    inline Vec3& operator-=(const Vec3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    // Negation
    inline Vec3 operator-() const {
        return Vec3(-x, -y, -z);
    }

    // Scalar multiplication
    inline Vec3 operator*(double scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }

    inline Vec3& operator*=(double scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    // Scalar division
    inline Vec3 operator/(double scalar) const {
        double inv = 1.0 / scalar;
        return Vec3(x * inv, y * inv, z * inv);
    }

    inline Vec3& operator/=(double scalar) {
        double inv = 1.0 / scalar;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }

    // Dot product
    inline double dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // Cross product
    inline Vec3 cross(const Vec3& other) const {
        return Vec3(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }

    // Squared length (avoid sqrt when possible)
    inline double length_squared() const {
        return x * x + y * y + z * z;
    }

    // Length (magnitude)
    inline double length() const {
        return std::sqrt(length_squared());
    }

    // Normalize in place
    inline Vec3& normalize() {
        double len = length();
        if (len > 1e-12) {
            double inv = 1.0 / len;
            x *= inv;
            y *= inv;
            z *= inv;
        }
        return *this;
    }

    // Return normalized copy
    inline Vec3 normalized() const {
        Vec3 result = *this;
        result.normalize();
        return result;
    }
};

// Free function: scalar * vec
inline Vec3 operator*(double scalar, const Vec3& v) {
    return Vec3(v.x * scalar, v.y * scalar, v.z * scalar);
}

// Free function: dot product
inline double dot(const Vec3& a, const Vec3& b) {
    return a.dot(b);
}

// Free function: cross product
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return a.cross(b);
}

// Free function: distance squared (avoid sqrt)
inline double distance_squared(const Vec3& a, const Vec3& b) {
    Vec3 diff = a - b;
    return diff.length_squared();
}

// Free function: distance
inline double distance(const Vec3& a, const Vec3& b) {
    return std::sqrt(distance_squared(a, b));
}

#endif // VEC3_HPP

