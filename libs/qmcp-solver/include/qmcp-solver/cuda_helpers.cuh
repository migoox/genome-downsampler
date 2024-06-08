#ifndef QMCP_SOLVER_CUDA_HELPERS
#define QMCP_SOLVER_CUDA_HELPERS

#include <iostream>

#include "cuda_runtime.h"

#define CUDA_HELPER_STRINGIFY(x) #x
// do while statement ensures that macro behaves like a single statement
#define CHECK_CUDA_ERROR(cuda_call)                                                   \
    do {                                                                              \
        cudaError_t cuda_status = cuda_call;                                          \
        if (cuda_status != cudaSuccess) {                                             \
            std::cerr << "CUDA error in " << CUDA_HELPER_STRINGIFY(cuda_call) << ": " \
                      << cudaGetErrorString(cuda_status) << std::endl;                \
            std::terminate();                                                         \
        }                                                                             \
    } while (0)

namespace qmcp {
namespace cuda {

template <class T>
inline T* malloc(size_t count) {
    T* ptr = nullptr;
    CHECK_CUDA_ERROR(cudaMalloc(reinterpret_cast<void**>(&ptr), count * sizeof(T)));
    return ptr;
}

template <class T>
inline void memcpy(T* dst, T* src, size_t count, cudaMemcpyKind kind) {
    CHECK_CUDA_ERROR(cudaMemcpy(dst, src, count * sizeof(T), kind));
}

// Copy from host to device
template <class T>
inline void memcpy_host_dev(T* dst, T* src, size_t count) {
    memcpy<T>(dst, src, count, cudaMemcpyHostToDevice);
}

// Copy from device to host
template <class T>
inline void memcpy_dev_host(T* dst, T* src, size_t count) {
    memcpy<T>(dst, src, count, cudaMemcpyDeviceToHost);
}

inline void free(void* dev_ptr) { CHECK_CUDA_ERROR(cudaFree(dev_ptr)); }

}  // namespace cuda
}  // namespace qmcp

#endif