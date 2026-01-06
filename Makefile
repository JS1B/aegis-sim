# Aegis-Murmuration Makefile
# Stage 1: Serial Foundation

CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic
OPTFLAGS := -O3 -march=native -ffast-math
DEBUGFLAGS := -g -fsanitize=address -DDEBUG

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
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

clean:
	rm -rf $(BUILDDIR) $(TARGET)

# Dependencies (auto-generated would be better, but keeping it simple)
$(BUILDDIR)/main.o: $(SRCDIR)/main.cpp $(INCDIR)/particle_system.hpp $(INCDIR)/flocking.hpp $(INCDIR)/vtk_writer.hpp
$(BUILDDIR)/particle_system.o: $(SRCDIR)/particle_system.cpp $(INCDIR)/particle_system.hpp $(INCDIR)/vec3.hpp $(INCDIR)/memory.hpp
$(BUILDDIR)/flocking.o: $(SRCDIR)/flocking.cpp $(INCDIR)/flocking.hpp $(INCDIR)/particle_system.hpp $(INCDIR)/vec3.hpp
$(BUILDDIR)/vtk_writer.o: $(SRCDIR)/vtk_writer.cpp $(INCDIR)/vtk_writer.hpp $(INCDIR)/particle_system.hpp

