#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include "flocking.hpp"
#include "particle_system.hpp"
#include "spatial_grid.hpp"
#include "vtk_writer.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

void print_usage(const char *program_name) {
  std::cerr << "Usage: " << program_name << " <N_particles> <timesteps> <T>\n";
  std::cerr << "\n";
  std::cerr << "Arguments:\n";
  std::cerr << "  N_particles     Number of particles to simulate\n";
  std::cerr << "  timesteps       Total number of simulation timesteps\n";
  std::cerr << "  T               Simulation timestep (dt, e.g., 0.01)\n";
  std::cerr << "\n";
  std::cerr << "Example:\n";
  std::cerr << "  " << program_name << " 50000 1000 0.01\n";
}

int main(int argc, char *argv[]) {
  // Parse command line arguments
  if (argc != 4) {
    print_usage(argv[0]);
    return 1;
  }

  size_t n_particles = static_cast<size_t>(std::atol(argv[1]));
  int timesteps = std::atoi(argv[2]);
  double dt = std::atof(argv[3]);

  if (n_particles == 0 || timesteps <= 0 || dt <= 0.0) {
    std::cerr << "Error: Invalid parameters\n";
    print_usage(argv[0]);
    return 1;
  }

  std::cout << "=== Aegis-Murmuration Simulation ===\n";
  std::cout << "Particles:       " << n_particles << "\n";
  std::cout << "Timesteps:       " << timesteps << "\n";
  std::cout << "Timestep (dt):   " << dt << "\n";

#ifdef _OPENMP
  int num_threads = omp_get_max_threads();
  std::cout << "OpenMP threads:  " << num_threads << "\n";
#else
  std::cout << "OpenMP:          disabled\n";
#endif

  std::cout << "\n";

  // Simulation parameters
  const int num_species = 3;    // Number of species for coloring
  const unsigned int seed = 42; // Random seed

  // Initialize particle system
  std::cout << "Initializing particle system..." << std::flush;
  auto init_start = std::chrono::high_resolution_clock::now();

  ParticleSystem particles(n_particles, num_species, seed);

  auto init_end = std::chrono::high_resolution_clock::now();
  double init_time =
      std::chrono::duration<double>(init_end - init_start).count();
  std::cout << " done (" << init_time << " s)\n";

  // Setup flocking parameters
  FlockingParams params;
  // Tune for spherical surface behavior - optimized for visible flocking
  params.separation_radius = 0.08;
  params.separation_weight = 2.5; // Strong separation to prevent collapse
  params.alignment_radius = 0.15;
  params.alignment_weight = 1.0; // Increased for coordinated movement
  params.cohesion_radius = 0.20;
  params.cohesion_weight = 1.2; // Increased to encourage grouping
  params.max_force = 5.0;       // Increased for more responsive behavior
  params.damping = 0.005;       // Reduced damping for more sustained movement

  // Allocate force buffer
  std::vector<Vec3> forces(n_particles);

  // Create spatial grid for optimized neighbor search
  // Grid resolution: 20 cells per dimension (tunable)
  const int grid_resolution = 20;
  const int grid_rebuild_interval = 10; // Rebuild every N steps

  std::cout << "Creating spatial grid (" << grid_resolution << "^3 cells)...\n";
  SpatialGrid grid(grid_resolution);

  // Build initial grid
  auto grid_start = std::chrono::high_resolution_clock::now();
  grid.build(particles);
  auto grid_end = std::chrono::high_resolution_clock::now();
  double grid_time =
      std::chrono::duration<double, std::milli>(grid_end - grid_start).count();

  auto stats = grid.get_stats();
  std::cout << "Grid built in " << grid_time << " ms\n";
  std::cout << "  Occupied cells: " << stats.occupied_cells << " / "
            << stats.total_cells << " ("
            << (100.0 * stats.occupied_cells / stats.total_cells) << "%)\n";
  std::cout << "  Avg particles/cell: " << std::fixed << std::setprecision(1)
            << stats.avg_particles_per_cell << "\n";
  std::cout << "  Max particles/cell: " << stats.max_particles_per_cell
            << "\n\n";

  // Write initial state
  std::cout << "Writing initial VTK mesh...\n";
  write_vtk_timestep(particles, "output_data/particles", 0, true);

  // Main simulation loop
  std::cout << "\nStarting simulation...\n";
  auto sim_start = std::chrono::high_resolution_clock::now();

  double total_grid_rebuild_time = 0.0;
  int grid_rebuild_count = 0;

  for (int t = 1; t <= timesteps; ++t) {
    // Rebuild grid periodically
    if (t % grid_rebuild_interval == 0) {
      auto rebuild_start = std::chrono::high_resolution_clock::now();
      grid.build(particles);
      auto rebuild_end = std::chrono::high_resolution_clock::now();
      total_grid_rebuild_time +=
          std::chrono::duration<double, std::milli>(rebuild_end - rebuild_start)
              .count();
      ++grid_rebuild_count;
    }

    // Compute flocking forces using spatial grid optimization
    compute_flocking_forces_grid(particles, params, grid, forces.data());

    // Integrate (applies spherical constraint) - now parallelized
    particles.integrate(forces.data(), dt);

    // Write VTK output for this timestep
    write_vtk_timestep(particles, "output_data/particles", t, false);

    // Output progress for each step
    auto now = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(now - sim_start).count();
    double steps_per_sec = (elapsed > 0) ? (t / elapsed) : 0.0;
    double eta = (steps_per_sec > 0) ? ((timesteps - t) / steps_per_sec) : 0.0;

    std::cout << "\rStep " << t << "/" << timesteps << " (" << std::fixed
              << std::setprecision(1) << (100.0 * t / timesteps) << "%) "
              << "| " << std::setprecision(0) << steps_per_sec << " steps/s "
              << "| ETA: " << std::setprecision(1) << eta << "s     "
              << std::flush;
  }

  auto sim_end = std::chrono::high_resolution_clock::now();
  double sim_time = std::chrono::duration<double>(sim_end - sim_start).count();

  // Compute final velocity statistics
  double min_vel = 1e10, max_vel = 0.0, avg_vel = 0.0;
  int moving_count = 0;
  for (size_t i = 0; i < n_particles; ++i) {
    Vec3 vel = particles.velocity(i);
    double vel_mag = vel.length();
    avg_vel += vel_mag;
    if (vel_mag < min_vel)
      min_vel = vel_mag;
    if (vel_mag > max_vel)
      max_vel = vel_mag;
    if (vel_mag > 0.001)
      moving_count++;
  }
  avg_vel /= n_particles;

  std::cout << "\n\n=== Simulation Complete ===\n";
  std::cout << "Total time:      " << std::fixed << std::setprecision(2)
            << sim_time << " s\n";
  std::cout << "Steps/second:    " << std::setprecision(1)
            << (timesteps / sim_time) << "\n";
  std::cout << "VTK files:       " << (timesteps + 1)
            << " (particles_000000.vtp to particles_" << std::setw(6)
            << std::setfill('0') << timesteps << ".vtp)\n";

  std::cout << "\n--- Final Velocity Statistics ---\n";
  std::cout << "Min velocity:    " << std::scientific << std::setprecision(3)
            << min_vel << "\n";
  std::cout << "Max velocity:    " << max_vel << "\n";
  std::cout << "Avg velocity:    " << avg_vel << "\n";
  std::cout << "Moving particles:" << std::fixed << std::setprecision(1)
            << moving_count << " (" << (100.0 * moving_count / n_particles)
            << "%)\n";

  if (grid_rebuild_count > 0) {
    std::cout << "\n--- Grid Statistics ---\n";
    std::cout << "Grid rebuilds:   " << grid_rebuild_count << "\n";
    std::cout << "Avg rebuild time:" << std::setprecision(2)
              << (total_grid_rebuild_time / grid_rebuild_count) << " ms\n";
    std::cout << "Total grid time: " << total_grid_rebuild_time << " ms\n";
  }

  std::cout << "\nOutput files written to: output_data/particles_*.vtp\n";
  std::cout << "Open in ParaView for visualization.\n";

  return 0;
}
