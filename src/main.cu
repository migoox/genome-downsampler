#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include <iostream>

const int kN = 10;  // Or any suitable size

// Kernel definition
__global__ void MatAdd(float A[kN][kN], float B[kN][kN], float C[kN][kN]) {
    int i = static_cast<int>(threadIdx.x);
    int j = static_cast<int>(threadIdx.y);
    C[i][j] = A[i][j] + B[i][j];
}

int main() {
    float a[kN][kN];
    float b[kN][kN];
    float c[kN][kN];

    float variable_two = 0.5F;  // NOLINT

    // Initialize matrices A and B with some values
    for (int i = 0; i < kN; ++i) {
        for (int j = 0; j < kN; ++j) {
            a[i][j] = static_cast<float>(i * kN + j);
            b[i][j] = static_cast<float>((i * kN + j) * 2);
        }
    }

    // Kernel invocation with one block of N * N * 1 threads
    int num_blocks = 1;
    dim3 threads_per_block(kN, kN);
    MatAdd<<<num_blocks, threads_per_block>>>(a, b, c);

    // Copy result from device memory to host memory
    cudaMemcpy(c, c, kN * kN * sizeof(float), cudaMemcpyDeviceToHost);

    // Display the result matrix C
    std::cout << "Result Matrix C:" << std::endl;
    for (auto& i : c) {
        for (float j : i) {
            std::cout << j << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
