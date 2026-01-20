#include "flocking.hpp"
#include "spatial_grid.hpp"
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

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

/**
 * Compute flocking force from a pre-computed neighbor list.
 * This is the optimized version used by the grid-based algorithm.
 */
static Vec3 compute_force_from_neighbors(
    const ParticleSystem& particles,
    const FlockingParams& params,
    size_t particle_index,
    const std::vector<size_t>& neighbors
) {
    Vec3 force;
    Vec3 pos_i = particles.position(particle_index);
    Vec3 vel_i = particles.velocity(particle_index);
    
    const double* px = particles.pos_x();
    const double* py = particles.pos_y();
    const double* pz = particles.pos_z();
    const double* vx = particles.vel_x();
    const double* vy = particles.vel_y();
    const double* vz = particles.vel_z();
    
    // Precompute squared radii
    double sep_radius_sq = params.separation_radius * params.separation_radius;
    double align_radius_sq = params.alignment_radius * params.alignment_radius;
    double coh_radius_sq = params.cohesion_radius * params.cohesion_radius;
    
    // Accumulators for each force component
    Vec3 separation_force;
    Vec3 avg_velocity;
    Vec3 center_of_mass;
    int align_count = 0;
    int coh_count = 0;
    
    // Single pass through neighbors for all three forces
    for (size_t j : neighbors) {
        if (j == particle_index) continue;
        
        // Compute distance
        double dx = pos_i.x - px[j];
        double dy = pos_i.y - py[j];
        double dz = pos_i.z - pz[j];
        double dist_sq = dx * dx + dy * dy + dz * dz;
        
        // Separation: strong repulsion at close range
        if (dist_sq < sep_radius_sq && dist_sq > 1e-12) {
            double dist = std::sqrt(dist_sq);
            double strength = params.separation_weight * (params.separation_radius - dist) / dist;
            separation_force.x += dx * strength;
            separation_force.y += dy * strength;
            separation_force.z += dz * strength;
        }
        
        // Alignment: match neighbor velocities
        if (dist_sq < align_radius_sq) {
            avg_velocity.x += vx[j];
            avg_velocity.y += vy[j];
            avg_velocity.z += vz[j];
            ++align_count;
        }
        
        // Cohesion: move toward local center of mass
        if (dist_sq < coh_radius_sq) {
            center_of_mass.x += px[j];
            center_of_mass.y += py[j];
            center_of_mass.z += pz[j];
            ++coh_count;
        }
    }
    
    // Add separation force
    force += separation_force;
    
    // Add alignment force
    if (align_count > 0) {
        avg_velocity /= static_cast<double>(align_count);
        force += (avg_velocity - vel_i) * params.alignment_weight;
    }
    
    // Add cohesion force
    if (coh_count > 0) {
        center_of_mass /= static_cast<double>(coh_count);
        force += (center_of_mass - pos_i) * params.cohesion_weight;
    }
    
    // Apply velocity damping
    force -= vel_i * params.damping;
    
    // Clamp force magnitude
    double force_mag_sq = force.length_squared();
    if (force_mag_sq > params.max_force * params.max_force) {
        double scale = params.max_force / std::sqrt(force_mag_sq);
        force *= scale;
    }
    
    return force;
}

void compute_flocking_forces_grid(
    const ParticleSystem& particles,
    const FlockingParams& params,
    const SpatialGrid& grid,
    Vec3* forces
) {
    size_t n = particles.count();
    
    // Find maximum radius for neighbor queries
    double max_radius = std::max({
        params.separation_radius,
        params.alignment_radius,
        params.cohesion_radius
    });
    
    // Parallel force computation using grid-based neighbor search
    #pragma omp parallel
    {
        // Thread-local neighbor buffer to avoid allocations
        std::vector<size_t> neighbors;
        neighbors.reserve(200);  // Pre-allocate for typical neighbor count
        
        #pragma omp for schedule(guided)
        for (size_t i = 0; i < n; ++i) {
            Vec3 pos_i = particles.position(i);
            
            // Query neighbors using spatial grid
            grid.query_neighbors(pos_i, max_radius, particles, neighbors);
            
            // Compute force from neighbor list
            forces[i] = compute_force_from_neighbors(particles, params, i, neighbors);
        }
    }
}

