# Aegis-Murmuration

WIP

**High-Performance Spherical Swarm & Dynamic Voronoi Tessellation**

A simulation of hundreds of thousands of agents "dancing" on a sphere, forming a dynamic honeycomb (Voronoi) shield. This project demonstrates HPC techniques including cache-optimized data structures, spherical geometry, and flocking algorithms.

## Stage 2: Spatial Partitioning & OpenMP Parallelization ✅

This implementation includes:

- **Structure of Arrays (SoA)** particle system for cache-efficient access
- **64-byte aligned memory** for AVX-512 compatibility
- **Spherical constraint** - particles are projected onto the unit sphere
- **3D Spatial Grid** - Voxel-based neighbor search (20³ cells)
- **OpenMP parallelization** - Force computation, integration, and grid building
- **Grid-based flocking** - O(N × k) instead of O(N²) with Separation, Alignment, and Cohesion
- **Adaptive grid rebuild** - Rebuild every 10 timesteps for optimal performance
- **Spherical Delaunay triangulation** via Qhull convex hull
- **Binary VTK mesh export** for ParaView visualization with connected triangles
- **Performance**: 12-260 steps/s for 50K-10K particles on 24 cores (see `STAGE2_PERFORMANCE.md`)

## How the Mesh Works

For points constrained to a unit sphere, the **3D convex hull IS the spherical Delaunay triangulation**:

- All points lie on the sphere surface
- Convex hull only connects surface points
- No edges or triangles pass through the sphere interior
- Result: Geographic tessellation like Earth's tectonic plates

This satisfies the **Euler characteristic**: V - E + F = 2

- For N vertices: expect ~2N triangles and ~3N edges

## Dependencies

- **Qhull** library for convex hull computation

```bash
# Arch Linux
sudo pacman -S qhull

# Ubuntu/Debian  
sudo apt install libqhull-dev

# HPC cluster
module load qhull/2020.2
```

## Building

```bash
# Release build (optimized with OpenMP)
make release

# Debug build (with sanitizers)
make debug

# Clean build artifacts
make clean
```

## Running

### Local Execution

```bash
# Default parameters: 5000 particles, 500 steps, output every 10 steps
./run_app.sh

# Custom parameters
./run_app.sh <N_particles> <timesteps> <output_interval>
./run_app.sh 50000 1000 10
```

### HPC Cluster (DTU/LSF)

```bash
# Submit job to HPC queue
bsub < run_app_hpc.sh
```

## Output

VTK PolyData mesh files are written to `output_data/`.

Run `python visualization/render_animation.py` to render an animation of the simulation. The frames are written to `visualization/frames/`.

## Performance Notes

- **Flocking forces**: O(N × k) with spatial grid, ~12 steps/s for 50K particles (24 cores)
- **Grid construction**: O(N) parallel build, ~14ms for 50K particles
- **Mesh building**: O(N log N) via Qhull, ~150ms for 10K particles
- **Point data**: O(N²) but OpenMP parallel, ~70ms for 10K particles
- **OpenMP scaling**: 6× speedup on 24 cores (10K particles)

**See `STAGE2_PERFORMANCE.md` for detailed benchmarks.**

## Future Stages

- **Stage 3**: GPU acceleration with CUDA/OpenCL + Jump Flooding Algorithm
- **Stage 4**: Advanced ParaView visualization + Voronoi cell coloring

## License

MIT
