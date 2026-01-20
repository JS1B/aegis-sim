#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include "particle_system.hpp"
#include "mesh.hpp"
#include <string>

/**
 * Write particle system as spherical Delaunay mesh to VTK PolyData (.vtp).
 * 
 * Computes convex hull (spherical Delaunay) and exports connected triangles.
 * Uses binary encoding with base64 for compact output.
 * 
 * The output includes:
 * - Point coordinates
 * - Triangle connectivity (Polys)
 * - Velocity vectors (for glyph visualization)
 * - Velocity magnitude (for scalar coloring)
 * - Local density (for coloring by neighbor count)
 * - Species ID (for categorical coloring)
 * 
 * @param particles The particle system to export
 * @param filename Output filename (should end in .vtp)
 * @param verbose Print mesh statistics
 * @return true on success, false on failure
 */
bool write_vtk(
    const ParticleSystem& particles,
    const std::string& filename,
    bool verbose = false
);

/**
 * Write precomputed mesh to VTK PolyData (.vtp).
 * Use this when you want to reuse the same mesh across multiple exports.
 * 
 * @param particles The particle system (for coordinates and velocities)
 * @param mesh Precomputed mesh with triangles and point data
 * @param filename Output filename
 * @return true on success, false on failure
 */
bool write_vtk_mesh(
    const ParticleSystem& particles,
    const Mesh& mesh,
    const std::string& filename
);

/**
 * Write particle system to VTK PolyData with timestep in filename.
 * Creates files like: prefix_000123.vtp
 * 
 * @param particles The particle system to export
 * @param prefix Filename prefix (e.g., "output_data/particles")
 * @param timestep Current timestep number
 * @param verbose Print mesh statistics
 * @return true on success, false on failure
 */
bool write_vtk_timestep(
    const ParticleSystem& particles,
    const std::string& prefix,
    int timestep,
    bool verbose = false
);

#endif // VTK_WRITER_HPP
