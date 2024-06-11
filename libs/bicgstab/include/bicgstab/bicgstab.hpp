#include <ctype.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <typeinfo>  // for usage of C++ typeid

#include "cublas_v2.h"
#include "cusparse_v2.h"
#include "helper_cuda.h"
#include "helper_cusolver.h"
#include "mmio.h"
#include "mmio_wrapper.h"

// profiling the code
#define TIME_INDIVIDUAL_LIBRARY_CALLS

#define DBICGSTAB_MAX_ULP_ERR 100
#define DBICGSTAB_EPS 1.E-14f  // 9e-2

#define CLEANUP()                                                             \
    do {                                                                      \
        if (x) free(x);                                                       \
        if (f) free(f);                                                       \
        if (r) free(r);                                                       \
        if (rw) free(rw);                                                     \
        if (p) free(p);                                                       \
        if (pw) free(pw);                                                     \
        if (s) free(s);                                                       \
        if (t) free(t);                                                       \
        if (v) free(v);                                                       \
        if (tx) free(tx);                                                     \
        if (Aval) free(Aval);                                                 \
        if (AcolsIndex) free(AcolsIndex);                                     \
        if (ArowsIndex) free(ArowsIndex);                                     \
        if (Mval) free(Mval);                                                 \
        if (devPtrX) checkCudaErrors(cudaFree(devPtrX));                      \
        if (devPtrF) checkCudaErrors(cudaFree(devPtrF));                      \
        if (devPtrR) checkCudaErrors(cudaFree(devPtrR));                      \
        if (devPtrRW) checkCudaErrors(cudaFree(devPtrRW));                    \
        if (devPtrP) checkCudaErrors(cudaFree(devPtrP));                      \
        if (devPtrS) checkCudaErrors(cudaFree(devPtrS));                      \
        if (devPtrT) checkCudaErrors(cudaFree(devPtrT));                      \
        if (devPtrV) checkCudaErrors(cudaFree(devPtrV));                      \
        if (devPtrAval) checkCudaErrors(cudaFree(devPtrAval));                \
        if (devPtrAcolsIndex) checkCudaErrors(cudaFree(devPtrAcolsIndex));    \
        if (devPtrArowsIndex) checkCudaErrors(cudaFree(devPtrArowsIndex));    \
        if (devPtrMval) checkCudaErrors(cudaFree(devPtrMval));                \
        if (stream) checkCudaErrors(cudaStreamDestroy(stream));               \
        if (cublasHandle) checkCudaErrors(cublasDestroy(cublasHandle));       \
        if (cusparseHandle) checkCudaErrors(cusparseDestroy(cusparseHandle)); \
        fflush(stdout);                                                       \
    } while (0)

static void gpu_pbicgstab(
    cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle, int m, int n, int nnz,
    const cusparseMatDescr_t descra, /* the coefficient matrix in CSR format */
    double* a, int* ia, int* ja,
    const cusparseMatDescr_t
        descrm, /* the preconditioner in CSR format, lower & upper triangular factor */
    double* vm, int* im, int* jm, cusparseSolveAnalysisInfo_t info_l,
    cusparseSolveAnalysisInfo_t info_u, /* the analysis of the lower and upper triangular parts */
    double* f, double* r, double* rw, double* p, double* pw, double* s, double* t, double* v,
    double* x, int maxit, double tol, double ttt_sv);