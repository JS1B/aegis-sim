#include "mesh_builder.hpp"
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullVertex.h>
#include <stdexcept>
#include <iostream>
#include <chrono>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

void compute_velocity_magnitudes(
    const ParticleSystem& particles,
    std::vector<double>& magnitudes
) {
    size_t n = particles.count();
    magnitudes.resize(n);
    
    const double* vx = particles.vel_x();
    const double* vy = particles.vel_y();
    const double* vz = particles.vel_z();
    
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; ++i) {
        magnitudes[i] = std::sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
    }
}

void compute_local_densities(
    const ParticleSystem& particles,
    double radius,
    std::vector<double>& densities
) {
    size_t n = particles.count();
    densities.resize(n);
    
    const double* px = particles.pos_x();
    const double* py = particles.pos_y();
    const double* pz = particles.pos_z();
    
    double radius_sq = radius * radius;
    
    // Brute-force O(N^2) but parallelized
    #pragma omp parallel for schedule(guided)
    for (size_t i = 0; i < n; ++i) {
        int count = 0;
        double xi = px[i], yi = py[i], zi = pz[i];
        
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            
            double dx = xi - px[j];
            double dy = yi - py[j];
            double dz = zi - pz[j];
            double dist_sq = dx * dx + dy * dy + dz * dz;
            
            if (dist_sq < radius_sq) {
                ++count;
            }
        }
        
        densities[i] = static_cast<double>(count);
    }
}

Mesh build_spherical_mesh(
    const ParticleSystem& particles,
    const MeshBuilderParams& params
) {
    size_t n = particles.count();
    
    if (n < 4) {
        throw std::runtime_error("Need at least 4 points for convex hull");
    }
    
    Mesh mesh;
    mesh.num_vertices = n;
    mesh.reserve(n);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Prepare point data for Qhull (interleaved x,y,z format)
    std::vector<double> points(n * 3);
    const double* px = particles.pos_x();
    const double* py = particles.pos_y();
    const double* pz = particles.pos_z();
    
    for (size_t i = 0; i < n; ++i) {
        points[i * 3 + 0] = px[i];
        points[i * 3 + 1] = py[i];
        points[i * 3 + 2] = pz[i];
    }
    
    // Compute convex hull using Qhull
    // Options: Qt = triangulate output (ensure all facets are triangles)
    orgQhull::Qhull qhull;
    
    try {
        qhull.runQhull("", 3, static_cast<int>(n), points.data(), "Qt");
    } catch (const orgQhull::QhullError& e) {
        throw std::runtime_error(std::string("Qhull error: ") + e.what());
    }
    
    // Extract triangular facets
    for (const orgQhull::QhullFacet& facet : qhull.facetList()) {
        if (facet.isGood()) {
            orgQhull::QhullVertexSet vertices = facet.vertices();
            
            if (vertices.size() == 3) {
                // Get vertex indices (Qhull point IDs)
                auto it = vertices.begin();
                int32_t v0 = (*it).point().id(); ++it;
                int32_t v1 = (*it).point().id(); ++it;
                int32_t v2 = (*it).point().id();
                
                mesh.add_triangle(v0, v1, v2);
            }
        }
    }
    
    auto hull_time = std::chrono::high_resolution_clock::now();
    
    // Compute point data in parallel
    compute_velocity_magnitudes(particles, mesh.velocity_magnitudes);
    compute_local_densities(particles, params.density_radius, mesh.local_densities);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    
    if (params.verbose) {
        double hull_ms = std::chrono::duration<double, std::milli>(hull_time - start_time).count();
        double data_ms = std::chrono::duration<double, std::milli>(end_time - hull_time).count();
        
        std::cout << "Mesh built: " << mesh.num_triangles << " triangles from " 
                  << mesh.num_vertices << " vertices\n";
        std::cout << "  Hull time: " << hull_ms << " ms\n";
        std::cout << "  Data time: " << data_ms << " ms\n";
        std::cout << "  Edges (computed): " << mesh.compute_num_edges() << "\n";
        
        // Verify Euler characteristic
        if (mesh.verify_euler()) {
            std::cout << "  Euler check: PASSED (V-E+F=2)\n";
        } else {
            int64_t V = mesh.num_vertices;
            int64_t F = mesh.num_triangles;
            int64_t E = mesh.compute_num_edges();
            std::cout << "  Euler check: WARNING (V-E+F=" << (V - E + F) 
                      << ", expected 2)\n";
        }
    }
    
    return mesh;
}

