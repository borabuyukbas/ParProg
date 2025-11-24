# High-Performance N-Body Simulation

A high-performance astronomical simulation capable of simulating over **1 million bodies** in real-time. This project leverages **CUDA** for GPU acceleration and **OpenMP** for multi-core CPU processing, implementing the efficient **Barnes-Hut algorithm** ($O(n \log n)$) to handle massive scale interactions.

## Key Features

*   **Massive Scale**: Simulates 1,000,000+ astronomical objects.
*   **High Performance**:
    *   **GPU Acceleration**: Custom CUDA kernels for massive parallelism.
    *   **CPU Parallelism**: OpenMP for multi-threaded execution.
    *   **Algorithmic Optimization**: Barnes-Hut tree algorithm for efficient force calculation.
*   **Visualization**: Generates high-resolution video visualizations of the universe evolution.
*   **Modern C++**: Built with C++20 standards.

## Technical Stack

*   **Languages**: C++20, CUDA C++
*   **Parallel Computing**: NVIDIA CUDA, OpenMP
*   **Build System**: CMake
*   **Tools**: FFmpeg (visualization), CLI11

## Quick Start

### Prerequisites
*   CMake (3.10+)
*   C++20 Compiler
*   CUDA Toolkit (11.0+)

### Build & Run
The project is organized into labs. `lab_3` contains the final, most advanced implementation.

```bash
cd lab_3
mkdir build && cd build
cmake .. && make

# Run simulation with CUDA (Mode 4) for 1000 bodies
./source/lab --num-bodies 1000 --simulation-mode 4
```

### Simulation Modes
*   `0`: Sequential (Baseline)
*   `1`: Parallel CPU (OpenMP)
*   `2`: Barnes-Hut (Tree-based)
*   `4`: **CUDA (GPU Accelerated)**

## Visualization
The simulation outputs frame data which can be rendered into video:
```bash
# Inside the build directory after running a simulation
../create_mp4.sh
```
*Generates `out.mp4` showing the universe evolution.*
