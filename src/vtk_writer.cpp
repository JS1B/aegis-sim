#include "vtk_writer.hpp"
#include "mesh_builder.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <vector>
#include <cstdint>

// Base64 encoding table
static const char base64_chars[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789+/";

/**
 * Encode binary data to base64.
 */
static std::string base64_encode(const unsigned char* data, size_t len) {
    std::string result;
    result.reserve(((len + 2) / 3) * 4);
    
    for (size_t i = 0; i < len; i += 3) {
        uint32_t n = static_cast<uint32_t>(data[i]) << 16;
        if (i + 1 < len) n |= static_cast<uint32_t>(data[i + 1]) << 8;
        if (i + 2 < len) n |= static_cast<uint32_t>(data[i + 2]);
        
        result.push_back(base64_chars[(n >> 18) & 0x3F]);
        result.push_back(base64_chars[(n >> 12) & 0x3F]);
        result.push_back((i + 1 < len) ? base64_chars[(n >> 6) & 0x3F] : '=');
        result.push_back((i + 2 < len) ? base64_chars[n & 0x3F] : '=');
    }
    
    return result;
}

/**
 * Encode array with header (VTK format: 4-byte size header + data).
 */
template<typename T>
static std::string encode_array_base64(const T* data, size_t count) {
    // VTK binary format: header (uint32 size in bytes) + raw data
    uint32_t byte_size = static_cast<uint32_t>(count * sizeof(T));
    
    std::vector<unsigned char> buffer(sizeof(uint32_t) + byte_size);
    
    // Copy header
    std::memcpy(buffer.data(), &byte_size, sizeof(uint32_t));
    
    // Copy data
    std::memcpy(buffer.data() + sizeof(uint32_t), data, byte_size);
    
    return base64_encode(buffer.data(), buffer.size());
}

bool write_vtk_mesh(
    const ParticleSystem& particles,
    const Mesh& mesh,
    const std::string& filename
) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    size_t n = particles.count();
    size_t num_polys = mesh.num_triangles;
    
    // Prepare interleaved position data (x0,y0,z0, x1,y1,z1, ...)
    std::vector<double> positions(n * 3);
    for (size_t i = 0; i < n; ++i) {
        positions[i * 3 + 0] = particles.pos_x()[i];
        positions[i * 3 + 1] = particles.pos_y()[i];
        positions[i * 3 + 2] = particles.pos_z()[i];
    }
    
    // Prepare velocity vectors
    std::vector<double> velocities(n * 3);
    for (size_t i = 0; i < n; ++i) {
        velocities[i * 3 + 0] = particles.vel_x()[i];
        velocities[i * 3 + 1] = particles.vel_y()[i];
        velocities[i * 3 + 2] = particles.vel_z()[i];
    }
    
    // Prepare species data
    std::vector<int32_t> species(n);
    for (size_t i = 0; i < n; ++i) {
        species[i] = particles.species_id()[i];
    }
    
    // Write VTK XML PolyData format with Polys (triangles)
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
    file << "  <PolyData>\n";
    file << "    <Piece NumberOfPoints=\"" << n 
         << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"" 
         << num_polys << "\">\n";
    
    // Point data (attributes per vertex)
    file << "      <PointData Vectors=\"Velocity\" Scalars=\"VelocityMagnitude\">\n";
    
    // Velocity vectors
    file << "        <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"binary\">\n";
    file << "          " << encode_array_base64(velocities.data(), velocities.size()) << "\n";
    file << "        </DataArray>\n";
    
    // Velocity magnitude
    file << "        <DataArray type=\"Float64\" Name=\"VelocityMagnitude\" format=\"binary\">\n";
    file << "          " << encode_array_base64(mesh.velocity_magnitudes.data(), mesh.velocity_magnitudes.size()) << "\n";
    file << "        </DataArray>\n";
    
    // Local density
    file << "        <DataArray type=\"Float64\" Name=\"LocalDensity\" format=\"binary\">\n";
    file << "          " << encode_array_base64(mesh.local_densities.data(), mesh.local_densities.size()) << "\n";
    file << "        </DataArray>\n";
    
    // Species ID
    file << "        <DataArray type=\"Int32\" Name=\"SpeciesID\" format=\"binary\">\n";
    file << "          " << encode_array_base64(species.data(), species.size()) << "\n";
    file << "        </DataArray>\n";
    
    file << "      </PointData>\n";
    
    // Points (coordinates)
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">\n";
    file << "          " << encode_array_base64(positions.data(), positions.size()) << "\n";
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    
    // Polys (triangles from spherical Delaunay)
    file << "      <Polys>\n";
    
    // Connectivity: triangle vertex indices
    // Format: t0_v0, t0_v1, t0_v2, t1_v0, t1_v1, t1_v2, ...
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n";
    file << "          " << encode_array_base64(mesh.triangle_indices.data(), mesh.triangle_indices.size()) << "\n";
    file << "        </DataArray>\n";
    
    // Offsets: cumulative vertex count per polygon
    // For triangles: 3, 6, 9, 12, ...
    std::vector<int32_t> offsets(num_polys);
    for (size_t i = 0; i < num_polys; ++i) {
        offsets[i] = static_cast<int32_t>((i + 1) * 3);
    }
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n";
    file << "          " << encode_array_base64(offsets.data(), offsets.size()) << "\n";
    file << "        </DataArray>\n";
    
    file << "      </Polys>\n";
    
    file << "    </Piece>\n";
    file << "  </PolyData>\n";
    file << "</VTKFile>\n";
    
    return file.good();
}

bool write_vtk(
    const ParticleSystem& particles,
    const std::string& filename,
    bool verbose
) {
    // Build spherical Delaunay mesh
    MeshBuilderParams params;
    params.verbose = verbose;
    params.density_radius = 0.3;
    
    Mesh mesh = build_spherical_mesh(particles, params);
    
    return write_vtk_mesh(particles, mesh, filename);
}

bool write_vtk_timestep(
    const ParticleSystem& particles,
    const std::string& prefix,
    int timestep,
    bool verbose
) {
    std::ostringstream filename;
    filename << prefix << "_" << std::setw(6) << std::setfill('0') << timestep << ".vtp";
    return write_vtk(particles, filename.str(), verbose);
}
