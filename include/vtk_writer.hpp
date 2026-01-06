#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include "particle_system.hpp"
#include <string>

/**
 * Write particle system to VTK PolyData (.vtp) format.
 * Uses binary encoding with base64 for compact output.
 * 
 * The output includes:
 * - Point coordinates
 * - Velocity vectors (for glyph visualization)
 * - Species ID (for coloring)
 * - Velocity magnitude (for scalar coloring)
 * 
 * @param particles The particle system to export
 * @param filename Output filename (should end in .vtp)
 * @return true on success, false on failure
 */
bool write_vtk(const ParticleSystem& particles, const std::string& filename);

/**
 * Write particle system to VTK PolyData with timestep in filename.
 * Creates files like: prefix_000123.vtp
 * 
 * @param particles The particle system to export
 * @param prefix Filename prefix (e.g., "output_data/particles")
 * @param timestep Current timestep number
 * @return true on success, false on failure
 */
bool write_vtk_timestep(
    const ParticleSystem& particles,
    const std::string& prefix,
    int timestep
);

#endif // VTK_WRITER_HPP

