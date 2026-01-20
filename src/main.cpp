#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <vector>

#include "particle_system.hpp"
#include "flocking.hpp"
#include "vtk_writer.hpp"

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <N_particles> <timesteps> <output_interval>\n";
    std::cerr << "\n";
    std::cerr << "Arguments:\n";
    std::cerr << "  N_particles     Number of particles to simulate\n";
    std::cerr << "  timesteps       Total number of simulation timesteps\n";
    std::cerr << "  output_interval Write VTK output every N steps (0 = no output)\n";
    std::cerr << "\n";
    std::cerr << "Example:\n";
    std::cerr << "  " << program_name << " 50000 1000 10\n";
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    if (argc != 4) {
        print_usage(argv[0]);
        return 1;
    }
    
    size_t n_particles = static_cast<size_t>(std::atol(argv[1]));
    int timesteps = std::atoi(argv[2]);
    int output_interval = std::atoi(argv[3]);
    
    if (n_particles == 0 || timesteps <= 0) {
        std::cerr << "Error: Invalid parameters\n";
        print_usage(argv[0]);
        return 1;
    }
    
    std::cout << "=== Aegis-Murmuration Simulation ===\n";
    std::cout << "Particles:       " << n_particles << "\n";
    std::cout << "Timesteps:       " << timesteps << "\n";
    std::cout << "Output interval: " << output_interval << "\n";
    std::cout << "\n";
    
    // Simulation parameters
    const double dt = 0.01;  // Timestep
    const int num_species = 3;  // Number of species for coloring
    const unsigned int seed = 42;  // Random seed
    
    // Initialize particle system
    std::cout << "Initializing particle system..." << std::flush;
    auto init_start = std::chrono::high_resolution_clock::now();
    
    ParticleSystem particles(n_particles, num_species, seed);
    
    auto init_end = std::chrono::high_resolution_clock::now();
    double init_time = std::chrono::duration<double>(init_end - init_start).count();
    std::cout << " done (" << init_time << " s)\n";
    
    // Setup flocking parameters
    FlockingParams params;
    // Tune for spherical surface behavior
    params.separation_radius = 0.08;
    params.separation_weight = 3.0;  // Strong separation to prevent collapse
    params.alignment_radius = 0.15;
    params.alignment_weight = 0.3;
    params.cohesion_radius = 0.25;
    params.cohesion_weight = 0.5;
    params.max_force = 2.0;
    params.damping = 0.02;
    
    // Allocate force buffer
    std::vector<Vec3> forces(n_particles);
    
    // Write initial state (verbose to show mesh statistics)
    if (output_interval > 0) {
        std::cout << "Writing initial VTK mesh...\n";
        write_vtk_timestep(particles, "output_data/particles", 0, true);
    }
    
    // Main simulation loop
    std::cout << "\nStarting simulation...\n";
    auto sim_start = std::chrono::high_resolution_clock::now();
    
    int output_count = 1;
    for (int t = 1; t <= timesteps; ++t) {
        // Compute flocking forces
        compute_flocking_forces(particles, params, forces.data());
        
        // Integrate (applies spherical constraint)
        particles.integrate(forces.data(), dt);
        
        // Output progress
        if (t % 100 == 0 || t == timesteps) {
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(now - sim_start).count();
            double steps_per_sec = t / elapsed;
            double eta = (timesteps - t) / steps_per_sec;
            
            std::cout << "\rStep " << t << "/" << timesteps 
                      << " (" << std::fixed << std::setprecision(1) 
                      << (100.0 * t / timesteps) << "%) "
                      << "| " << std::setprecision(0) << steps_per_sec << " steps/s "
                      << "| ETA: " << std::setprecision(1) << eta << "s     " 
                      << std::flush;
        }
        
        // Write VTK output
        if (output_interval > 0 && t % output_interval == 0) {
            write_vtk_timestep(particles, "output_data/particles", output_count);
            ++output_count;
        }
    }
    
    auto sim_end = std::chrono::high_resolution_clock::now();
    double sim_time = std::chrono::duration<double>(sim_end - sim_start).count();
    
    std::cout << "\n\n=== Simulation Complete ===\n";
    std::cout << "Total time:      " << std::fixed << std::setprecision(2) << sim_time << " s\n";
    std::cout << "Steps/second:    " << std::setprecision(1) << (timesteps / sim_time) << "\n";
    std::cout << "VTK files:       " << output_count << "\n";
    
    if (output_interval > 0) {
        std::cout << "\nOutput files written to: output_data/particles_*.vtp\n";
        std::cout << "Open in ParaView for visualization.\n";
    }
    
    return 0;
}

