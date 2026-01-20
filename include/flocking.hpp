#ifndef FLOCKING_HPP
#define FLOCKING_HPP

#include "particle_system.hpp"
#include "vec3.hpp"

// Forward declaration to avoid circular dependency
class SpatialGrid;

/**
 * Flocking force parameters.
 * Tune these to control swarm behavior.
 */
struct FlockingParams {
    // Separation: Strong repulsion at close range
    double separation_radius = 0.05;
    double separation_weight = 2.0;
    
    // Alignment: Match neighbor velocities
    double alignment_radius = 0.2;
    double alignment_weight = 0.5;
    
    // Cohesion: Move toward local center of mass
    double cohesion_radius = 0.3;
    double cohesion_weight = 0.8;
    
    // Maximum force magnitude (prevents explosion)
    double max_force = 1.0;
    
    // Velocity damping (0 = no damping, 1 = full damping)
    double damping = 0.01;
};

/**
 * Compute flocking forces for all particles.
 * Uses brute-force O(N^2) neighbor search.
 * 
 * @param particles The particle system
 * @param params Flocking behavior parameters
 * @param forces Output array for computed forces (must be pre-allocated)
 */
void compute_flocking_forces(
    const ParticleSystem& particles,
    const FlockingParams& params,
    Vec3* forces
);

/**
 * Compute flocking force for a single particle.
 * Useful for debugging or when updating specific particles.
 * 
 * @param particles The particle system
 * @param params Flocking behavior parameters
 * @param particle_index Index of particle to compute force for
 * @return The computed force vector
 */
Vec3 compute_flocking_force_single(
    const ParticleSystem& particles,
    const FlockingParams& params,
    size_t particle_index
);

/**
 * Compute flocking forces using spatial grid optimization.
 * Uses O(N * avg_neighbors) instead of O(N^2) with grid-based neighbor search.
 * Parallelized with OpenMP for maximum performance.
 * 
 * @param particles The particle system
 * @param params Flocking behavior parameters
 * @param grid Spatial grid for neighbor queries
 * @param forces Output array for computed forces (must be pre-allocated)
 */
void compute_flocking_forces_grid(
    const ParticleSystem& particles,
    const FlockingParams& params,
    const SpatialGrid& grid,
    Vec3* forces
);

#endif // FLOCKING_HPP

