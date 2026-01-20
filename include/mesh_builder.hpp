#ifndef MESH_BUILDER_HPP
#define MESH_BUILDER_HPP

#include "mesh.hpp"
#include "particle_system.hpp"

/**
 * Parameters for mesh building and point data computation.
 */
struct MeshBuilderParams {
    // Radius for local density calculation
    double density_radius = 0.3;
    
    // Enable verbose output (Euler check, timing)
    bool verbose = false;
};

/**
 * Build spherical Delaunay mesh from particle positions.
 * 
 * Uses Qhull to compute the 3D convex hull of points on the unit sphere.
 * For points constrained to a sphere, the convex hull IS the spherical
 * Delaunay triangulation - no internal faces, only surface triangles.
 * 
 * This ensures no visual artifacts (edges through the sphere) in ParaView.
 * 
 * @param particles The particle system with positions/velocities
 * @param params Optional parameters for mesh building
 * @return Mesh with triangles and computed point data
 * @throws std::runtime_error if Qhull fails
 */
Mesh build_spherical_mesh(
    const ParticleSystem& particles,
    const MeshBuilderParams& params = MeshBuilderParams{}
);

/**
 * Compute velocity magnitudes for all particles.
 * Uses OpenMP for parallel computation.
 * 
 * @param particles The particle system
 * @param magnitudes Output array (must be pre-sized to particles.count())
 */
void compute_velocity_magnitudes(
    const ParticleSystem& particles,
    std::vector<double>& magnitudes
);

/**
 * Compute local density (neighbor count) for all particles.
 * Uses brute-force O(N^2) with OpenMP parallelization.
 * 
 * @param particles The particle system
 * @param radius Neighborhood radius for density calculation
 * @param densities Output array (must be pre-sized to particles.count())
 */
void compute_local_densities(
    const ParticleSystem& particles,
    double radius,
    std::vector<double>& densities
);

#endif // MESH_BUILDER_HPP

