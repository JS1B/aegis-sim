# Aegis-Murmuration Makefile
# Stage 1.5: Serial Foundation + Spherical Mesh Connectivity

CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic
OPTFLAGS := -O3 -march=native -ffast-math -fopenmp
DEBUGFLAGS := -g -fsanitize=address -DDEBUG

# Qhull libraries for convex hull / spherical Delaunay
LDFLAGS := -lqhullcpp -lqhull_r -fopenmp

# Directories
SRCDIR := src
INCDIR := include
BUILDDIR := build
TARGET := aegis_sim

# Source files
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(BUILDDIR)/%.o)

# Default: Release build
.PHONY: all clean debug release

all: release

release: CXXFLAGS += $(OPTFLAGS)
release: $(TARGET)

debug: CXXFLAGS += $(DEBUGFLAGS)
debug: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

clean:
	rm -rf $(BUILDDIR) $(TARGET)

# Dependencies
$(BUILDDIR)/main.o: $(SRCDIR)/main.cpp $(INCDIR)/particle_system.hpp $(INCDIR)/flocking.hpp $(INCDIR)/vtk_writer.hpp
$(BUILDDIR)/particle_system.o: $(SRCDIR)/particle_system.cpp $(INCDIR)/particle_system.hpp $(INCDIR)/vec3.hpp $(INCDIR)/memory.hpp
$(BUILDDIR)/flocking.o: $(SRCDIR)/flocking.cpp $(INCDIR)/flocking.hpp $(INCDIR)/particle_system.hpp $(INCDIR)/vec3.hpp
$(BUILDDIR)/vtk_writer.o: $(SRCDIR)/vtk_writer.cpp $(INCDIR)/vtk_writer.hpp $(INCDIR)/particle_system.hpp $(INCDIR)/mesh.hpp
$(BUILDDIR)/mesh_builder.o: $(SRCDIR)/mesh_builder.cpp $(INCDIR)/mesh_builder.hpp $(INCDIR)/mesh.hpp $(INCDIR)/particle_system.hpp
