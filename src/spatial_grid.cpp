#include "spatial_grid.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

SpatialGrid::SpatialGrid(int resolution)
    : resolution_(resolution),
      grid_min_(-1.2), // Slightly larger than unit sphere radius
      grid_max_(1.2) {

  if (resolution_ <= 0) {
    resolution_ = 20; // Default fallback
  }

  cell_size_ = (grid_max_ - grid_min_) / static_cast<double>(resolution_);

  // Allocate cells (resolution^3 total cells)
  size_t total_cells = static_cast<size_t>(resolution_) *
                       static_cast<size_t>(resolution_) *
                       static_cast<size_t>(resolution_);
  cells_.resize(total_cells);
}

void SpatialGrid::clear() {
  for (auto &cell : cells_) {
    cell.clear();
  }
}

bool SpatialGrid::position_to_cell_coords(const Vec3 &pos, int &x, int &y,
                                          int &z) const {
  // Map position to grid coordinates
  double fx = (pos.x - grid_min_) / cell_size_;
  double fy = (pos.y - grid_min_) / cell_size_;
  double fz = (pos.z - grid_min_) / cell_size_;

  x = static_cast<int>(std::floor(fx));
  y = static_cast<int>(std::floor(fy));
  z = static_cast<int>(std::floor(fz));

  // Check bounds
  if (x < 0 || x >= resolution_ || y < 0 || y >= resolution_ || z < 0 ||
      z >= resolution_) {
    return false;
  }

  return true;
}

int SpatialGrid::position_to_cell_index(const Vec3 &pos) const {
  int x, y, z;
  if (!position_to_cell_coords(pos, x, y, z)) {
    return -1;
  }
  return cell_index_3d_to_1d(x, y, z);
}

void SpatialGrid::build(const ParticleSystem &particles) {
  size_t n = particles.count();

  // Clear existing grid
  clear();

  // Two-phase parallel construction to avoid race conditions

  // Phase 1: Count particles per cell (parallel with atomics)
  std::vector<size_t> cell_counts(cells_.size(), 0);

#pragma omp parallel
  {
    // Thread-local counts to reduce atomic contention
    std::vector<size_t> local_counts(cells_.size(), 0);

#pragma omp for schedule(static) nowait
    for (size_t i = 0; i < n; ++i) {
      Vec3 pos = particles.position(i);
      int cell_idx = position_to_cell_index(pos);

      if (cell_idx >= 0) {
        local_counts[cell_idx]++;
      }
    }

// Reduce thread-local counts to global counts
#pragma omp critical
    {
      for (size_t i = 0; i < cells_.size(); ++i) {
        cell_counts[i] += local_counts[i];
      }
    }
  }

  // Phase 2: Reserve space in each cell
  for (size_t i = 0; i < cells_.size(); ++i) {
    if (cell_counts[i] > 0) {
      cells_[i].reserve(cell_counts[i]);
    }
  }

// Phase 3: Fill cells (parallel, using push_back which is now safe)
// We use thread-local vectors and merge to avoid synchronization
#ifdef _OPENMP
  int num_threads = omp_get_max_threads();
#else
  int num_threads = 1;
#endif

  // Create thread-local cell vectors
  std::vector<std::vector<std::vector<size_t>>> thread_cells(num_threads);
  for (int t = 0; t < num_threads; ++t) {
    thread_cells[t].resize(cells_.size());
  }

#pragma omp parallel
  {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif

    auto &local_cells = thread_cells[thread_id];

#pragma omp for schedule(static) nowait
    for (size_t i = 0; i < n; ++i) {
      Vec3 pos = particles.position(i);
      int cell_idx = position_to_cell_index(pos);

      if (cell_idx >= 0) {
        local_cells[cell_idx].push_back(i);
      }
    }
  }

// Merge thread-local cells into main grid
#pragma omp parallel for schedule(guided)
  for (size_t cell_idx = 0; cell_idx < cells_.size(); ++cell_idx) {
    for (int t = 0; t < num_threads; ++t) {
      auto &local_cell = thread_cells[t][cell_idx];
      cells_[cell_idx].insert(cells_[cell_idx].end(), local_cell.begin(),
                              local_cell.end());
    }
  }
}

void SpatialGrid::query_neighbors(const Vec3 &pos, double radius,
                                  const ParticleSystem &particles,
                                  std::vector<size_t> &neighbors) const {
  neighbors.clear();

  // Get cell coordinates for query position
  int cx, cy, cz;
  if (!position_to_cell_coords(pos, cx, cy, cz)) {
    // Position outside grid - fall back to checking all particles (rare)
    neighbors.reserve(particles.count());
    double radius_sq = radius * radius;

    for (size_t i = 0; i < particles.count(); ++i) {
      Vec3 p = particles.position(i);
      double dist_sq = (p.x - pos.x) * (p.x - pos.x) +
                       (p.y - pos.y) * (p.y - pos.y) +
                       (p.z - pos.z) * (p.z - pos.z);
      if (dist_sq < radius_sq) {
        neighbors.push_back(i);
      }
    }
    return;
  }

  // Calculate cell range to check (3x3x3 stencil, possibly larger for big
  // radius)
  int cell_radius = static_cast<int>(std::ceil(radius / cell_size_));

  double radius_sq = radius * radius;

  // Check all cells within range
  for (int dz = -cell_radius; dz <= cell_radius; ++dz) {
    int z = cx + dz;
    if (z < 0 || z >= resolution_)
      continue;

    for (int dy = -cell_radius; dy <= cell_radius; ++dy) {
      int y = cy + dy;
      if (y < 0 || y >= resolution_)
        continue;

      for (int dx = -cell_radius; dx <= cell_radius; ++dx) {
        int x = cx + dx;
        if (x < 0 || x >= resolution_)
          continue;

        int cell_idx = cell_index_3d_to_1d(x, y, z);
        const auto &cell = cells_[cell_idx];

        // Check each particle in this cell
        for (size_t particle_idx : cell) {
          Vec3 p = particles.position(particle_idx);
          double dist_sq = (p.x - pos.x) * (p.x - pos.x) +
                           (p.y - pos.y) * (p.y - pos.y) +
                           (p.z - pos.z) * (p.z - pos.z);

          if (dist_sq < radius_sq) {
            neighbors.push_back(particle_idx);
          }
        }
      }
    }
  }
}

SpatialGrid::GridStats SpatialGrid::get_stats() const {
  GridStats stats;
  stats.total_cells = cells_.size();
  stats.occupied_cells = 0;
  stats.total_particles = 0;
  stats.max_particles_per_cell = 0;
  stats.min_particles_per_cell = std::numeric_limits<size_t>::max();

  for (const auto &cell : cells_) {
    size_t count = cell.size();
    if (count > 0) {
      stats.occupied_cells++;
      stats.total_particles += count;
      stats.max_particles_per_cell =
          std::max(stats.max_particles_per_cell, count);
      stats.min_particles_per_cell =
          std::min(stats.min_particles_per_cell, count);
    }
  }

  if (stats.occupied_cells > 0) {
    stats.avg_particles_per_cell = static_cast<double>(stats.total_particles) /
                                   static_cast<double>(stats.occupied_cells);
  } else {
    stats.avg_particles_per_cell = 0.0;
    stats.min_particles_per_cell = 0;
  }

  return stats;
}
