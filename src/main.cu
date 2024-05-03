#include <htslib/hts.h>
#include <stdio.h>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "qmcp-solver/qmcp-solver.hpp"
#include "qmcp-solver/sequential_cost_scaling_network_solver.hpp"
cudaError_t addWithCuda(int* c, const int* a, const int* b, unsigned int size);

__global__ void addKernel(int* c, const int* a, const int* b) {
    int i = static_cast<int>(threadIdx.x);
    c[i] = a[i] + b[i];
}

int main() {
    // Bam api and qmcp solver test
    int M = 500;
    auto bam_path = std::filesystem::path(
        "/home/borys/Downloads/gpu-programming/data/ESIB_EQA_2023.SARS2.01/"
        "reads.bam");
    auto solver = qmcp::SequentialCostScalingNetworkSolver(M, bam_path);

    auto start = std::chrono::high_resolution_clock::now();
    solver.solve();
    auto stop = std::chrono::high_resolution_clock::now();
    auto solve_duration =
        std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    std::cout << "SOLVE TOOK " << solve_duration.count() << "[seconds]"
              << std::endl;

    // Define some variables
    const int array_size = 5;
    const int a[array_size] = {1, 2, 3, 4, 5};
    const int b[array_size] = {10, 20, 30, 40, 50};
    int c[array_size] = {0};

    // Add vectors in parallel.
    cudaError_t cuda_status = addWithCuda(c, a, b, array_size);
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");  // NOLINT
        return 1;
    }

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",  // NOLINT
           c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cuda_status = cudaDeviceReset();
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");  // NOLINT
        return 1;
    }

    return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int* c, const int* a, const int* b, unsigned int size) {
    int* dev_a = 0;
    int* dev_b = 0;
    int* dev_c = 0;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaError_t cuda_status = cudaSetDevice(0);
    if (cuda_status != cudaSuccess) {
        fprintf(  // NOLINT
            stderr,
            "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;  // NOLINT
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cuda_status = cudaMalloc((void**)&dev_c, size * sizeof(int));  // NOLINT
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");  // NOLINT
        goto Error;                             // NOLINT
    }

    cuda_status = cudaMalloc((void**)&dev_a, size * sizeof(int));  // NOLINT
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");  // NOLINT
        goto Error;                             // NOLINT
    }

    cuda_status = cudaMalloc((void**)&dev_b, size * sizeof(int));  // NOLINT
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");  // NOLINT
        goto Error;                             // NOLINT
    }

    // Copy input vectors from host memory to GPU buffers.
    cuda_status =
        cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");  // NOLINT
        goto Error;                             // NOLINT
    }

    cuda_status =
        cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cuda_status != cudaSuccess) {
        fprintf(stderr,  // NOLINT
                "cudaMemcpy failed!");
        goto Error;  // NOLINT
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cuda_status = cudaGetLastError();
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n",  // NOLINT
                cudaGetErrorString(cuda_status));
        goto Error;  // NOLINT
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cuda_status = cudaDeviceSynchronize();
    if (cuda_status != cudaSuccess) {
        fprintf(stderr,  // NOLINT
                "cudaDeviceSynchronize returned error code %d after launching "
                "addKernel!\n",
                cuda_status);
        goto Error;  // NOLINT
    }

    // Copy output vector from GPU buffer to host memory.
    cuda_status =
        cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cuda_status != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");  // NOLINT
        goto Error;                             // NOLINT
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);

    return cuda_status;
}
