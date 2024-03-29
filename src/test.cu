#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <iostream>

const int N = 10;// Or any suitable size

// Kernel definition
__global__ void MatAdd(float A[N][N], float B[N][N], float C[N][N])
{
    int i = threadIdx.x;
    int j = threadIdx.y;
    C[i][j] = A[i][j] + B[i][j];
}

int main()
{
    float A[N][N], B[N][N], C[N][N];

    // Initialize matrices A and B with some values
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = i * N + j;// Example initialization, you can use any method
            B[i][j] = (i * N + j) * 2;// Example initialization, you can use any method
        }
    }

    // Kernel invocation with one block of N * N * 1 threads
    int numBlocks = 1;
    dim3 threadsPerBlock(N, N);
    MatAdd<<<numBlocks, threadsPerBlock>>>(A, B, C);

    // Copy result from device memory to host memory
    cudaMemcpy(C, C, N * N * sizeof(float), cudaMemcpyDeviceToHost);

    // Display the result matrix C
    std::cout << "Result Matrix C:" << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) { std::cout << C[i][j] << " "; }
        std::cout << std::endl;
    }

    return 0;
}
