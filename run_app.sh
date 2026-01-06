#!/bin/bash
# Local test run for Aegis-Murmuration simulation

BINARY="./aegis_sim"

# Check if binary exists
if [ ! -f "$BINARY" ]; then
    echo "Binary not found! Building..."
    make release
    if [ $? -ne 0 ]; then
        echo "Build failed!"
        exit 1
    fi
fi

# Create output directory
mkdir -p output_data

# Default parameters if not provided
N_PARTICLES=${1:-5000}
TIMESTEPS=${2:-500}
OUTPUT_INTERVAL=${3:-10}

echo "Running simulation with:"
echo "  Particles: $N_PARTICLES"
echo "  Timesteps: $TIMESTEPS"
echo "  Output interval: $OUTPUT_INTERVAL"
echo ""

# Run simulation
$BINARY $N_PARTICLES $TIMESTEPS $OUTPUT_INTERVAL

