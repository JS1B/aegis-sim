#ifndef PARTICLE_SYSTEM_HPP
#define PARTICLE_SYSTEM_HPP

#include <cstddef>
#include "memory.hpp"
#include "vec3.hpp"

/**
 * Structure of Arrays (SoA) particle system for cache-efficient access.
 * All particles are constrained to the unit sphere.
 */
class ParticleSystem {
public:
    /**
     * Create a particle system with N particles.
     * Particles are randomly initialized on the unit sphere.
     * @param count Number of particles
     * @param num_species Number of species for coloring (default: 1)
     * @param seed Random seed for reproducibility
     */
    ParticleSystem(size_t count, int num_species = 1, unsigned int seed = 42);
    
    ~ParticleSystem();
    
    // No copy
    ParticleSystem(const ParticleSystem&) = delete;
    ParticleSystem& operator=(const ParticleSystem&) = delete;
    
    // Move
    ParticleSystem(ParticleSystem&& other) noexcept;
    ParticleSystem& operator=(ParticleSystem&& other) noexcept;
    
    /**
     * Get position of particle i as Vec3.
     */
    Vec3 position(size_t i) const;
    
    /**
     * Set position of particle i from Vec3.
     */
    void set_position(size_t i, const Vec3& pos);
    
    /**
     * Get velocity of particle i as Vec3.
     */
    Vec3 velocity(size_t i) const;
    
    /**
     * Set velocity of particle i from Vec3.
     */
    void set_velocity(size_t i, const Vec3& vel);
    
    /**
     * Apply spherical constraint: normalize position to unit sphere
     * and project velocity to be tangent to sphere.
     */
    void apply_spherical_constraint(size_t i);
    
    /**
     * Apply spherical constraint to all particles.
     */
    void apply_spherical_constraint_all();
    
    /**
     * Integrate velocities and positions with given timestep.
     * Applies spherical constraint after integration.
     * @param forces Array of force vectors (one per particle)
     * @param dt Timestep
     */
    void integrate(const Vec3* forces, double dt);
    
    // Accessors
    size_t count() const { return count_; }
    
    // Direct array access for performance-critical loops
    double* pos_x() { return pos_x_; }
    double* pos_y() { return pos_y_; }
    double* pos_z() { return pos_z_; }
    double* vel_x() { return vel_x_; }
    double* vel_y() { return vel_y_; }
    double* vel_z() { return vel_z_; }
    int* species_id() { return species_id_; }
    
    const double* pos_x() const { return pos_x_; }
    const double* pos_y() const { return pos_y_; }
    const double* pos_z() const { return pos_z_; }
    const double* vel_x() const { return vel_x_; }
    const double* vel_y() const { return vel_y_; }
    const double* vel_z() const { return vel_z_; }
    const int* species_id() const { return species_id_; }

private:
    size_t count_;
    
    // Positions (aligned for SIMD)
    double* pos_x_;
    double* pos_y_;
    double* pos_z_;
    
    // Velocities (aligned for SIMD)
    double* vel_x_;
    double* vel_y_;
    double* vel_z_;
    
    // Species ID for multi-species flocking
    int* species_id_;
    
    /**
     * Initialize particles randomly on unit sphere.
     */
    void initialize_random(int num_species, unsigned int seed);
    
    /**
     * Free all allocated memory.
     */
    void free_memory();
};

#endif // PARTICLE_SYSTEM_HPP

