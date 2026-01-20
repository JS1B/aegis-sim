#ifndef SPATIAL_GRID_HPP
#define SPATIAL_GRID_HPP

#include "particle_system.hpp"
#include "vec3.hpp"
#include <vector>
#include <cstddef>

/**
 * 3D spatial grid for efficient neighbor queries.
 * 
 * Divides the 3D space containing the unit sphere into uniform voxel cells.
 * Particles are assigned to cells based on their positions, enabling
 * O(avg_neighbors) neighbor queries instead of O(N) brute-force search.
 * 
 * Grid construction is parallelized using OpenMP.
 */
class SpatialGrid {
public:
    /**
     * Create a spatial grid with specified resolution.
     * 
     * @param resolution Number of cells per dimension (default: 20)
     *                  Total cells = resolution^3
     */
    explicit SpatialGrid(int resolution = 20);
    
    /**
     * Build/rebuild the grid from current particle positions.
     * Uses parallel construction with OpenMP.
     * 
     * @param particles The particle system to index
     */
    void build(const ParticleSystem& particles);
    
    /**
     * Query all particles within a given radius of a position.
     * Only checks relevant cells (3x3x3 stencil around query position).
     * 
     * @param pos Query position
     * @param radius Search radius
     * @param particles The particle system (for position lookups)
     * @param neighbors Output vector of particle indices (cleared and filled)
     */
    void query_neighbors(
        const Vec3& pos,
        double radius,
        const ParticleSystem& particles,
        std::vector<size_t>& neighbors
    ) const;
    
    /**
     * Clear all grid cells (call before rebuild).
     */
    void clear();
    
    /**
     * Get grid statistics for debugging/profiling.
     */
    struct GridStats {
        size_t total_cells;
        size_t occupied_cells;
        size_t total_particles;
        double avg_particles_per_cell;
        size_t max_particles_per_cell;
        size_t min_particles_per_cell;
    };
    
    GridStats get_stats() const;
    
    // Accessors
    int resolution() const { return resolution_; }
    double cell_size() const { return cell_size_; }
    size_t num_cells() const { return cells_.size(); }

private:
    int resolution_;        // Cells per dimension
    double cell_size_;      // Size of each cell
    double grid_min_;       // Minimum coordinate (e.g., -1.2)
    double grid_max_;       // Maximum coordinate (e.g., 1.2)
    
    // Particle indices stored in each cell
    // cells_[cell_index] = vector of particle indices in that cell
    std::vector<std::vector<size_t>> cells_;
    
    /**
     * Convert 3D cell coordinates to 1D array index.
     * Uses row-major ordering: index = z * res^2 + y * res + x
     */
    inline int cell_index_3d_to_1d(int x, int y, int z) const {
        return z * resolution_ * resolution_ + y * resolution_ + x;
    }
    
    /**
     * Convert 3D position to cell coordinates.
     * Returns false if position is outside grid bounds.
     */
    bool position_to_cell_coords(const Vec3& pos, int& x, int& y, int& z) const;
    
    /**
     * Convert 3D position to 1D cell index.
     * Returns -1 if position is outside grid bounds.
     */
    int position_to_cell_index(const Vec3& pos) const;
    
    /**
     * Clamp cell coordinate to valid range [0, resolution-1].
     */
    inline int clamp_cell_coord(int coord) const {
        if (coord < 0) return 0;
        if (coord >= resolution_) return resolution_ - 1;
        return coord;
    }
};

#endif // SPATIAL_GRID_HPP
