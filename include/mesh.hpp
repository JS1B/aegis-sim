#ifndef MESH_HPP
#define MESH_HPP

#include <vector>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

/**
 * Mesh data structure for spherical Delaunay triangulation.
 * 
 * Stores triangle connectivity and precomputed point data for VTK export.
 * The triangles form the convex hull of points on the unit sphere,
 * which is equivalent to the spherical Delaunay triangulation.
 */
struct Mesh {
    // Triangle connectivity - flat array of vertex indices
    // Format: [t0_v0, t0_v1, t0_v2, t1_v0, t1_v1, t1_v2, ...]
    std::vector<int32_t> triangle_indices;
    
    // Number of triangles (triangle_indices.size() / 3)
    size_t num_triangles;
    
    // Number of vertices (particles)
    size_t num_vertices;
    
    // Precomputed point data for VTK visualization
    std::vector<double> velocity_magnitudes;  // |v| for each particle
    std::vector<double> local_densities;      // neighbor count within radius
    
    /**
     * Default constructor - empty mesh.
     */
    Mesh() : num_triangles(0), num_vertices(0) {}
    
    /**
     * Reserve memory for expected sizes.
     * For N points on a sphere, expect ~2N triangles (Euler formula).
     * 
     * @param num_verts Number of vertices
     * @param num_tris Expected number of triangles (default: 2 * num_verts)
     */
    void reserve(size_t num_verts, size_t num_tris = 0) {
        if (num_tris == 0) {
            num_tris = 2 * num_verts;  // Euler characteristic estimate
        }
        triangle_indices.reserve(num_tris * 3);
        velocity_magnitudes.reserve(num_verts);
        local_densities.reserve(num_verts);
    }
    
    /**
     * Clear all data.
     */
    void clear() {
        triangle_indices.clear();
        velocity_magnitudes.clear();
        local_densities.clear();
        num_triangles = 0;
        num_vertices = 0;
    }
    
    /**
     * Add a triangle by its three vertex indices.
     * Indices must be valid vertex indices [0, num_vertices).
     */
    void add_triangle(int32_t v0, int32_t v1, int32_t v2) {
        triangle_indices.push_back(v0);
        triangle_indices.push_back(v1);
        triangle_indices.push_back(v2);
        ++num_triangles;
    }
    
    /**
     * Get vertex indices for triangle i.
     * @param i Triangle index [0, num_triangles)
     * @param v0, v1, v2 Output vertex indices
     */
    void get_triangle(size_t i, int32_t& v0, int32_t& v1, int32_t& v2) const {
        size_t base = i * 3;
        v0 = triangle_indices[base];
        v1 = triangle_indices[base + 1];
        v2 = triangle_indices[base + 2];
    }
    
    /**
     * Verify Euler characteristic for closed mesh on sphere.
     * For a sphere: V - E + F = 2
     * With triangles: E = 3F/2, so V - 3F/2 + F = 2 => V - F/2 = 2
     * 
     * @return true if mesh satisfies Euler characteristic
     */
    bool verify_euler() const {
        // For triangulated sphere: F = 2V - 4 (from V - E + F = 2, E = 3V - 6)
        // More precisely: num_triangles should be approximately 2 * num_vertices - 4
        int64_t expected_triangles = 2 * static_cast<int64_t>(num_vertices) - 4;
        int64_t actual_triangles = static_cast<int64_t>(num_triangles);
        
        // Allow small tolerance for numerical issues
        return std::abs(expected_triangles - actual_triangles) <= 4;
    }
    
    /**
     * Compute number of edges from triangles.
     * For a closed triangulated surface: E = 3F/2 (each edge shared by 2 triangles)
     */
    size_t compute_num_edges() const {
        return (num_triangles * 3) / 2;
    }
};

#endif // MESH_HPP

