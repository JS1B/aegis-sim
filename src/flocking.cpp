#include "flocking.hpp"
#include <cmath>
#include <algorithm>

/**
 * Compute separation force - strong repulsion at close range.
 * Uses inverse distance for smooth falloff.
 */
static Vec3 compute_separation(
    const ParticleSystem& particles,
    size_t i,
    double radius,
    double weight
) {
    Vec3 force;
    Vec3 pos_i = particles.position(i);
    
    const double* px = particles.pos_x();
    const double* py = particles.pos_y();
    const double* pz = particles.pos_z();
    
    double radius_sq = radius * radius;
    
    for (size_t j = 0; j < particles.count(); ++j) {
        if (i == j) continue;
        
        // Use squared distance to avoid sqrt
        double dx = pos_i.x - px[j];
        double dy = pos_i.y - py[j];
        double dz = pos_i.z - pz[j];
        double dist_sq = dx * dx + dy * dy + dz * dz;
        
        if (dist_sq < radius_sq && dist_sq > 1e-12) {
            double dist = std::sqrt(dist_sq);
            // Inverse distance weighting: closer = stronger repulsion
            double strength = weight * (radius - dist) / dist;
            force.x += dx * strength;
            force.y += dy * strength;
            force.z += dz * strength;
        }
    }
    
    return force;
}

/**
 * Compute alignment force - match neighbor velocities.
 */
static Vec3 compute_alignment(
    const ParticleSystem& particles,
    size_t i,
    double radius,
    double weight
) {
    Vec3 avg_velocity;
    Vec3 pos_i = particles.position(i);
    Vec3 vel_i = particles.velocity(i);
    
    const double* px = particles.pos_x();
    const double* py = particles.pos_y();
    const double* pz = particles.pos_z();
    const double* vx = particles.vel_x();
    const double* vy = particles.vel_y();
    const double* vz = particles.vel_z();
    
    double radius_sq = radius * radius;
    int neighbor_count = 0;
    
    for (size_t j = 0; j < particles.count(); ++j) {
        if (i == j) continue;
        
        double dx = pos_i.x - px[j];
        double dy = pos_i.y - py[j];
        double dz = pos_i.z - pz[j];
        double dist_sq = dx * dx + dy * dy + dz * dz;
        
        if (dist_sq < radius_sq) {
            avg_velocity.x += vx[j];
            avg_velocity.y += vy[j];
            avg_velocity.z += vz[j];
            ++neighbor_count;
        }
    }
    
    if (neighbor_count > 0) {
        avg_velocity /= static_cast<double>(neighbor_count);
        // Steer towards average velocity
        return (avg_velocity - vel_i) * weight;
    }
    
    return Vec3();
}

/**
 * Compute cohesion force - move toward local center of mass.
 */
static Vec3 compute_cohesion(
    const ParticleSystem& particles,
    size_t i,
    double radius,
    double weight
) {
    Vec3 center_of_mass;
    Vec3 pos_i = particles.position(i);
    
    const double* px = particles.pos_x();
    const double* py = particles.pos_y();
    const double* pz = particles.pos_z();
    
    double radius_sq = radius * radius;
    int neighbor_count = 0;
    
    for (size_t j = 0; j < particles.count(); ++j) {
        if (i == j) continue;
        
        double dx = pos_i.x - px[j];
        double dy = pos_i.y - py[j];
        double dz = pos_i.z - pz[j];
        double dist_sq = dx * dx + dy * dy + dz * dz;
        
        if (dist_sq < radius_sq) {
            center_of_mass.x += px[j];
            center_of_mass.y += py[j];
            center_of_mass.z += pz[j];
            ++neighbor_count;
        }
    }
    
    if (neighbor_count > 0) {
        center_of_mass /= static_cast<double>(neighbor_count);
        // Steer towards center of mass
        return (center_of_mass - pos_i) * weight;
    }
    
    return Vec3();
}

Vec3 compute_flocking_force_single(
    const ParticleSystem& particles,
    const FlockingParams& params,
    size_t particle_index
) {
    Vec3 force;
    
    // Separation - strongest force at close range
    force += compute_separation(
        particles, particle_index,
        params.separation_radius, params.separation_weight
    );
    
    // Alignment - match neighbor velocities
    force += compute_alignment(
        particles, particle_index,
        params.alignment_radius, params.alignment_weight
    );
    
    // Cohesion - move toward local center
    force += compute_cohesion(
        particles, particle_index,
        params.cohesion_radius, params.cohesion_weight
    );
    
    // Apply velocity damping
    Vec3 vel = particles.velocity(particle_index);
    force -= vel * params.damping;
    
    // Clamp force magnitude
    double force_mag_sq = force.length_squared();
    if (force_mag_sq > params.max_force * params.max_force) {
        double scale = params.max_force / std::sqrt(force_mag_sq);
        force *= scale;
    }
    
    return force;
}

void compute_flocking_forces(
    const ParticleSystem& particles,
    const FlockingParams& params,
    Vec3* forces
) {
    // Brute-force O(N^2): compute force for each particle
    for (size_t i = 0; i < particles.count(); ++i) {
        forces[i] = compute_flocking_force_single(particles, params, i);
    }
}

