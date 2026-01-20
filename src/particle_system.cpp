#include "particle_system.hpp"
#include <random>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

ParticleSystem::ParticleSystem(size_t count, int num_species, unsigned int seed)
    : count_(count),
      pos_x_(nullptr), pos_y_(nullptr), pos_z_(nullptr),
      vel_x_(nullptr), vel_y_(nullptr), vel_z_(nullptr),
      species_id_(nullptr) {
    
    // Allocate aligned memory
    pos_x_ = aligned_alloc_array<double>(count_);
    pos_y_ = aligned_alloc_array<double>(count_);
    pos_z_ = aligned_alloc_array<double>(count_);
    vel_x_ = aligned_alloc_array<double>(count_);
    vel_y_ = aligned_alloc_array<double>(count_);
    vel_z_ = aligned_alloc_array<double>(count_);
    species_id_ = aligned_alloc_array<int>(count_);
    
    initialize_random(num_species, seed);
}

ParticleSystem::~ParticleSystem() {
    free_memory();
}

ParticleSystem::ParticleSystem(ParticleSystem&& other) noexcept
    : count_(other.count_),
      pos_x_(other.pos_x_), pos_y_(other.pos_y_), pos_z_(other.pos_z_),
      vel_x_(other.vel_x_), vel_y_(other.vel_y_), vel_z_(other.vel_z_),
      species_id_(other.species_id_) {
    other.count_ = 0;
    other.pos_x_ = other.pos_y_ = other.pos_z_ = nullptr;
    other.vel_x_ = other.vel_y_ = other.vel_z_ = nullptr;
    other.species_id_ = nullptr;
}

ParticleSystem& ParticleSystem::operator=(ParticleSystem&& other) noexcept {
    if (this != &other) {
        free_memory();
        count_ = other.count_;
        pos_x_ = other.pos_x_;
        pos_y_ = other.pos_y_;
        pos_z_ = other.pos_z_;
        vel_x_ = other.vel_x_;
        vel_y_ = other.vel_y_;
        vel_z_ = other.vel_z_;
        species_id_ = other.species_id_;
        
        other.count_ = 0;
        other.pos_x_ = other.pos_y_ = other.pos_z_ = nullptr;
        other.vel_x_ = other.vel_y_ = other.vel_z_ = nullptr;
        other.species_id_ = nullptr;
    }
    return *this;
}

void ParticleSystem::free_memory() {
    if (pos_x_) aligned_free_array(pos_x_);
    if (pos_y_) aligned_free_array(pos_y_);
    if (pos_z_) aligned_free_array(pos_z_);
    if (vel_x_) aligned_free_array(vel_x_);
    if (vel_y_) aligned_free_array(vel_y_);
    if (vel_z_) aligned_free_array(vel_z_);
    if (species_id_) aligned_free_array(species_id_);
}

void ParticleSystem::initialize_random(int num_species, unsigned int seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist_theta(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> dist_cos_phi(-1.0, 1.0);
    std::uniform_real_distribution<double> dist_vel(-0.1, 0.1);  // Increased from ±0.01 to ±0.1
    std::uniform_int_distribution<int> dist_species(0, num_species - 1);
    
    for (size_t i = 0; i < count_; ++i) {
        // Random point on unit sphere using spherical coordinates
        double theta = dist_theta(rng);      // Azimuthal angle [0, 2π]
        double cos_phi = dist_cos_phi(rng);  // cos(polar angle) [-1, 1]
        double sin_phi = std::sqrt(1.0 - cos_phi * cos_phi);
        
        pos_x_[i] = sin_phi * std::cos(theta);
        pos_y_[i] = sin_phi * std::sin(theta);
        pos_z_[i] = cos_phi;
        
        // Small random initial velocities (tangent to sphere)
        Vec3 pos(pos_x_[i], pos_y_[i], pos_z_[i]);
        Vec3 random_vel(dist_vel(rng), dist_vel(rng), dist_vel(rng));
        
        // Project velocity to be tangent to sphere: v_tan = v - (v · n) * n
        double radial_component = random_vel.dot(pos);
        vel_x_[i] = random_vel.x - radial_component * pos.x;
        vel_y_[i] = random_vel.y - radial_component * pos.y;
        vel_z_[i] = random_vel.z - radial_component * pos.z;
        
        // Assign species
        species_id_[i] = dist_species(rng);
    }
}

Vec3 ParticleSystem::position(size_t i) const {
    return Vec3(pos_x_[i], pos_y_[i], pos_z_[i]);
}

void ParticleSystem::set_position(size_t i, const Vec3& pos) {
    pos_x_[i] = pos.x;
    pos_y_[i] = pos.y;
    pos_z_[i] = pos.z;
}

Vec3 ParticleSystem::velocity(size_t i) const {
    return Vec3(vel_x_[i], vel_y_[i], vel_z_[i]);
}

void ParticleSystem::set_velocity(size_t i, const Vec3& vel) {
    vel_x_[i] = vel.x;
    vel_y_[i] = vel.y;
    vel_z_[i] = vel.z;
}

void ParticleSystem::apply_spherical_constraint(size_t i) {
    // Normalize position to unit sphere
    Vec3 pos(pos_x_[i], pos_y_[i], pos_z_[i]);
    double len = pos.length();
    
    if (len > 1e-12) {
        double inv_len = 1.0 / len;
        pos_x_[i] *= inv_len;
        pos_y_[i] *= inv_len;
        pos_z_[i] *= inv_len;
        
        // Project velocity to be tangent to sphere
        // v_tan = v - (v · n) * n, where n is the normalized position
        Vec3 n(pos_x_[i], pos_y_[i], pos_z_[i]);
        Vec3 vel(vel_x_[i], vel_y_[i], vel_z_[i]);
        double radial_component = vel.dot(n);
        
        vel_x_[i] -= radial_component * n.x;
        vel_y_[i] -= radial_component * n.y;
        vel_z_[i] -= radial_component * n.z;
    }
}

void ParticleSystem::apply_spherical_constraint_all() {
    for (size_t i = 0; i < count_; ++i) {
        apply_spherical_constraint(i);
    }
}

void ParticleSystem::integrate(const Vec3* forces, double dt) {
    // Parallel integration - each particle update is independent
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < count_; ++i) {
        // Update velocity: v += F * dt
        vel_x_[i] += forces[i].x * dt;
        vel_y_[i] += forces[i].y * dt;
        vel_z_[i] += forces[i].z * dt;
        
        // Update position: p += v * dt
        pos_x_[i] += vel_x_[i] * dt;
        pos_y_[i] += vel_y_[i] * dt;
        pos_z_[i] += vel_z_[i] * dt;
        
        // Apply spherical constraint
        apply_spherical_constraint(i);
    }
}

