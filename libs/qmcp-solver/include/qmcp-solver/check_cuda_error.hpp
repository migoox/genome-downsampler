#define STRINGIFY(x) #x
// do while statement ensures that macro behaves like a single statement
#define CHECK_CUDA_ERROR(cuda_call)                                       \
    do {                                                                  \
        cudaError_t cuda_status = cuda_call;                              \
        if (cuda_status != cudaSuccess) {                                 \
            std::cerr << "CUDA error in " << STRINGIFY(cuda_call) << ": " \
                      << cudaGetErrorString(cuda_status) << std::endl;    \
            std::terminate();                                             \
        }                                                                 \
    } while (0)