#include "naive_cuda_simulation.cuh"
#include "physics/gravitation.h"
#include "physics/mechanics.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cuda_wrappers.cuh"
#include "constants.h"

std::vector<double2> vector_map_vector2d_to_double2(const std::vector<Vector2d<double>> vector) {
    std::vector<double2> returnVector (vector.size());

    for (size_t i = 0; i < vector.size(); ++i) {
        returnVector[i] = double2{vector[i][0], vector[i][1]};
    }

    return returnVector;
}

void NaiveCudaSimulation::allocate_device_memory(Universe& universe, void** d_weights, void** d_forces, void** d_velocities, void** d_positions){

    size_t number_bodies = universe.num_bodies;
    size_t memory_size_weights = number_bodies * sizeof(double);
    size_t memory_size_vectors = number_bodies * sizeof(double2);

    // alloc weight memory 
    parprog_cudaMalloc(d_weights, memory_size_weights);

    // alloc vector memory
    parprog_cudaMalloc(d_forces, memory_size_vectors);
    parprog_cudaMalloc(d_velocities, memory_size_vectors);
    parprog_cudaMalloc(d_positions, memory_size_vectors);
}

void NaiveCudaSimulation::free_device_memory(void** d_weights, void** d_forces, void** d_velocities, void** d_positions){

    // free memory 
    parprog_cudaFree(*d_weights);
    parprog_cudaFree(*d_forces);
    parprog_cudaFree(*d_velocities);
    parprog_cudaFree(*d_positions);

    // avoid dangling pointers 
    d_weights = nullptr;
    d_forces = nullptr;
    d_velocities = nullptr;
    d_positions = nullptr;
}

void NaiveCudaSimulation::copy_data_to_device(Universe& universe, void* d_weights, void* d_forces, void* d_velocities, void* d_positions){
    parprog_cudaMemcpy(d_weights, universe.weights.data(), universe.num_bodies * sizeof(double), cudaMemcpyHostToDevice);
    parprog_cudaMemcpy(d_forces, vector_map_vector2d_to_double2(universe.forces).data(), universe.num_bodies * sizeof(double2), cudaMemcpyHostToDevice);
    parprog_cudaMemcpy(d_velocities, vector_map_vector2d_to_double2(universe.velocities).data(), universe.num_bodies * sizeof(double2), cudaMemcpyHostToDevice);
    parprog_cudaMemcpy(d_positions, vector_map_vector2d_to_double2(universe.positions).data(), universe.num_bodies * sizeof(double2), cudaMemcpyHostToDevice);
}

void NaiveCudaSimulation::copy_data_from_device(Universe& universe, void* d_weights, void* d_forces, void* d_velocities, void* d_positions){
    parprog_cudaMemcpy(universe.weights.data(), d_weights, universe.num_bodies * sizeof(double), cudaMemcpyDeviceToHost);

    double2* forces = (double2*) malloc(universe.num_bodies * sizeof(double2));
    parprog_cudaMemcpy(forces, d_forces, universe.num_bodies * sizeof(double2), cudaMemcpyDeviceToHost);
    for (int i = 0; i < universe.num_bodies; ++i) {
        double2 cur = forces[i];
        universe.forces[i] = Vector2d<double>(cur.x, cur.y);
    }
    free(forces);

    double2* velocities = (double2*) malloc(universe.num_bodies * sizeof(double2));
    parprog_cudaMemcpy(velocities, d_velocities, universe.num_bodies * sizeof(double2), cudaMemcpyDeviceToHost);
    for (int i = 0; i < universe.num_bodies; ++i) {
        double2 cur = velocities[i];
        universe.velocities[i] = Vector2d<double>(cur.x, cur.y);
    }
    free(velocities);

    double2* positions = (double2*) malloc(universe.num_bodies * sizeof(double2));
    parprog_cudaMemcpy(positions, d_positions, universe.num_bodies * sizeof(double2), cudaMemcpyDeviceToHost);
    for (int i = 0; i < universe.num_bodies; ++i) {
        double2 cur = positions[i];
        universe.positions[i] = Vector2d<double>(cur.x, cur.y);
    }
    free(positions);
}

__global__
void calculate_forces_kernel(std::uint32_t num_bodies, double2* d_positions, double* d_weights, double2* d_forces){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_bodies) return;

    double2 i_position = d_positions[i];
    double i_weight = d_weights[i];

    double2 total_force {0, 0};
    for (int j = 0; j < num_bodies; ++j) {
        if (j == i) continue;

        double2 j_position = d_positions[j];
        double j_weight = d_weights[j];
        double2 direction {j_position.x - i_position.x, j_position.y - i_position.y};

        double distance = sqrt(pow(direction.x, 2) + pow(direction.y, 2));
        double force = gravitational_constant * ((i_weight * j_weight)/(pow(distance, 2)));
        double unit_vector_force = force / distance;

        total_force.x += direction.x * unit_vector_force;
        total_force.y += direction.y * unit_vector_force;
    }

    d_forces[i] = total_force;
}

void NaiveCudaSimulation::calculate_forces(Universe& universe, void* d_positions, void* d_weights, void* d_forces){
    int block_size = 512;
    int grid_size = universe.num_bodies % block_size == 0 ? universe.num_bodies / block_size : (universe.num_bodies - (universe.num_bodies % block_size) + block_size) / block_size;
    
    dim3 block_dim(block_size);
    dim3 grid_dim(grid_size);
    calculate_forces_kernel<<<grid_dim, block_dim>>>(universe.num_bodies, (double2*) d_positions, (double*) d_weights, (double2*) d_forces);;
}

__global__
void calculate_velocities_kernel(std::uint32_t num_bodies, double2* d_forces, double* d_weights, double2* d_velocities){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_bodies) return;

    double m = d_weights[i];
    double2 force = d_forces[i];
    double2 v0 = d_velocities[i];
    double2 a = {force.x / m, force.y / m};
    d_velocities[i] = {v0.x + a.x * epoch_in_seconds, v0.y + a.y * epoch_in_seconds};
}

void NaiveCudaSimulation::calculate_velocities(Universe& universe, void* d_forces, void* d_weights, void* d_velocities){
    int block_size = 512;
    int grid_size = universe.num_bodies % block_size == 0 ? universe.num_bodies / block_size : (universe.num_bodies - (universe.num_bodies % block_size) + block_size) / block_size;

    dim3 block_dim(block_size);
    dim3 grid_dim(grid_size);

    calculate_velocities_kernel<<<grid_dim, block_dim>>>(universe.num_bodies, (double2 *)d_forces, (double *)d_weights, (double2 *)d_velocities);
}

__global__
void calculate_positions_kernel(std::uint32_t num_bodies, double2* d_velocities, double2* d_positions){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_bodies) return;

    double2 pos = d_positions[i];
    double2 vel = d_velocities[i];
    d_positions[i] = {pos.x + vel.x * epoch_in_seconds, pos.y + vel.y * epoch_in_seconds};
}

void NaiveCudaSimulation::calculate_positions(Universe& universe, void* d_velocities, void* d_positions){
    int block_size = 512;
    int grid_size = universe.num_bodies % block_size == 0 ? universe.num_bodies / block_size : (universe.num_bodies - (universe.num_bodies % block_size) + block_size) / block_size;

    dim3 block_dim(block_size);
    dim3 grid_dim(grid_size);
    calculate_positions_kernel<<<grid_dim, block_dim>>>(universe.num_bodies, (double2 *)d_velocities, (double2 *)d_positions);
}

void NaiveCudaSimulation::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){

    void* d_weights;
    void* d_forces;
    void* d_velocities;
    void* d_positions;

    allocate_device_memory(universe, &d_weights, &d_forces, &d_velocities, &d_positions);
    
    for (int i = 0; i < num_epochs; i++)
    {
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs, d_weights, d_forces, d_velocities, d_positions);
    }

    free_device_memory(d_weights, d_forces, d_velocities, d_positions);
}

__global__
void get_pixels_kernel(std::uint32_t num_bodies, double2* d_positions, std::uint8_t* d_pixels, std::uint32_t plot_width, std::uint32_t plot_height, double plot_bounding_box_x_min, double plot_bounding_box_x_max, double plot_bounding_box_y_min, double plot_bounding_box_y_max){
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // check if within bounding box
    if (d_positions[idx].x >= plot_bounding_box_x_min && d_positions[idx].x <= plot_bounding_box_x_max && d_positions[idx].y >= plot_bounding_box_y_min && d_positions[idx].y <= plot_bounding_box_y_max)
    {
        int pixel_x = (d_positions[idx].x - plot_bounding_box_x_min) / (plot_bounding_box_x_max - plot_bounding_box_x_min) * (plot_width - 1);
        int pixel_y = (d_positions[idx].y - plot_bounding_box_y_min) / (plot_bounding_box_y_max - plot_bounding_box_y_min) * (plot_height - 1);

        // write to 255 to corresponding pixel

        d_pixels[pixel_y * plot_width + pixel_x] = 255;
    }
    // just keep the zero
}

std::vector<std::uint8_t> NaiveCudaSimulation::get_pixels(std::uint32_t plot_width, std::uint32_t plot_height, BoundingBox plot_bounding_box, void* d_positions, std::uint32_t num_bodies){
    // allocate memory
    void* d_pixels;
    uint32_t number_pixels = plot_width * plot_height;
    parprog_cudaMalloc(&d_pixels, number_pixels * sizeof(uint8_t));

    std::vector<std::uint8_t> pixels;
    pixels.resize(number_pixels, 0);

    // call get_pixels_kernel (write either 1 or zero)

    dim3 blockDim(num_bodies, 1, 1);
    dim3 gridDim(1, 1);
    get_pixels_kernel<<<gridDim, blockDim>>>(num_bodies, reinterpret_cast<double2*>(d_positions), reinterpret_cast<uint8_t*>(d_pixels),
        plot_width, plot_height, plot_bounding_box.x_min, plot_bounding_box.x_max, plot_bounding_box.y_min, plot_bounding_box.y_max);

    // copy back from device
    parprog_cudaMemcpy(pixels.data(), &d_pixels, plot_width * plot_height * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    // free memory
    parprog_cudaFree(d_pixels);

    return pixels;
}

__global__
void compress_pixels_kernel(std::uint32_t num_raw_pixels, std::uint8_t* d_raw_pixels, std::uint8_t* d_compressed_pixels){
    // basically reduction algorithm
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // for every compressed pixel
    uint8_t sum = 0;
    // iterate over 8 pixels
    for (uint8_t i = 0; i < 8; i++)
    {
        if (idx * 8 + i >= num_raw_pixels)
            break;
        if (d_raw_pixels[idx * 8 + i] != 0)
        {
            sum += (1 << i);
        }
    }

    d_compressed_pixels[idx] = sum;
}

void NaiveCudaSimulation::compress_pixels(std::vector<std::uint8_t>& raw_pixels, std::vector<std::uint8_t>& compressed_pixels){
    
    // allocate memory
    void* d_raw_pixels;
    void* d_compressed_pixels;

    size_t number_raw_pixels = raw_pixels.size();
    size_t number_comp_pixels = compressed_pixels.size();

    int block_size = 512;
    int grid_size = number_compressed_pixels % block_size == 0 ? number_compressed_pixels / block_size : (number_compressed_pixels - (number_compressed_pixels % block_size) + block_size) / block_size;

    parprog_cudaMalloc(&d_raw_pixels, number_raw_pixels * sizeof(uint8_t));
    parprog_cudaMalloc(&d_compressed_pixels, number_comp_pixels * sizeof(uint8_t));

    // copy to device
    parprog_cudaMemcpy(d_raw_pixels, raw_pixels.data(), number_raw_pixels * sizeof(uint8_t), cudaMemcpyHostToDevice);
    // parprog_cudaMemcpy(d_compressed_pixels, &raw_pixels, number_comp_pixels * sizeof(uint8_t), cudaMemcpyHostToDevice);

    dim3 blockDim(block_size);
    dim3 gridDim(grid_size);

    // call kernel
    compress_pixels_kernel<<<gridDim, blockDim>>>(number_raw_pixels, reinterpret_cast<uint8_t*>(d_raw_pixels), reinterpret_cast<uint8_t*>(d_compressed_pixels));

    // copy back to host
    parprog_cudaMemcpy(compressed_pixels.data(), d_compressed_pixels, number_comp_pixels * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    // free memory
    parprog_cudaFree(d_raw_pixels);
    parprog_cudaFree(d_compressed_pixels);
}

void NaiveCudaSimulation::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs, void* d_weights, void* d_forces, void* d_velocities, void* d_positions){
    calculate_forces(universe, d_positions, d_weights, d_forces);
    calculate_velocities(universe, d_forces, d_weights, d_velocities);
    calculate_positions(universe, d_velocities, d_positions);

    universe.current_simulation_epoch++;
    if(create_intermediate_plots){
        if(universe.current_simulation_epoch % plot_intermediate_epochs == 0){
            std::vector<std::uint8_t> pixels = get_pixels(plotter.get_plot_width(), plotter.get_plot_height(), plotter.get_plot_bounding_box(), d_positions, universe.num_bodies);
            plotter.add_active_pixels_to_image(pixels);

            // This is a dummy to use compression in plotting, although not beneficial performance-wise
            // ----
            // std::vector<std::uint8_t> compressed_pixels;
            // compressed_pixels.resize(pixels.size()/8);
            // compress_pixels(pixels, compressed_pixels);
            // plotter.add_compressed_pixels_to_image(compressed_pixels);
            // ----

            plotter.write_and_clear();
        }
    }
}

void NaiveCudaSimulation::calculate_forces_kernel_test_adapter(std::uint32_t grid_dim, std::uint32_t block_dim, std::uint32_t num_bodies, void* d_positions, void* d_weights, void* d_forces){
    // adapter function used by automatic tests. DO NOT MODIFY.
    dim3 blockDim(block_dim);
    dim3 gridDim(grid_dim);
    calculate_forces_kernel<<<gridDim, blockDim>>>(num_bodies, (double2*) d_positions, (double*) d_weights, (double2*) d_forces);
}

void NaiveCudaSimulation::calculate_velocities_kernel_test_adapter(std::uint32_t grid_dim, std::uint32_t block_dim, std::uint32_t num_bodies, void* d_forces, void* d_weights, void* d_velocities){
    // adapter function used by automatic tests. DO NOT MODIFY.
    dim3 blockDim(block_dim);
    dim3 gridDim(grid_dim);
    calculate_velocities_kernel<<<gridDim, blockDim>>>(num_bodies, (double2*) d_forces, (double*) d_weights, (double2*) d_velocities);
}

void NaiveCudaSimulation::calculate_positions_kernel_test_adapter(std::uint32_t grid_dim, std::uint32_t block_dim, std::uint32_t num_bodies, void* d_velocities, void* d_positions){
    // adapter function used by automatic tests. DO NOT MODIFY.
    dim3 blockDim(block_dim);
    dim3 gridDim(grid_dim);
    calculate_positions_kernel<<<gridDim, blockDim>>>(num_bodies, (double2*) d_velocities, (double2*) d_positions);
}
