/*
 * Copyright 1993-2022 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO LICENSEE:
 *
 * This source code and/or documentation ("Licensed Deliverables") are
 * subject to NVIDIA intellectual property rights under U.S. and
 * international Copyright laws.
 *
 * These Licensed Deliverables contained herein is PROPRIETARY and
 * CONFIDENTIAL to NVIDIA and is being provided under the terms and
 * conditions of a form of NVIDIA software license agreement by and
 * between NVIDIA and Licensee ("License Agreement") or electronically
 * accepted by Licensee.  Notwithstanding any terms or conditions to
 * the contrary in the License Agreement, reproduction or disclosure
 * of the Licensed Deliverables to any third party without the express
 * written consent of NVIDIA is prohibited.
 *
 * NOTWITHSTANDING ANY TERMS OR CONDITIONS TO THE CONTRARY IN THE
 * LICENSE AGREEMENT, NVIDIA MAKES NO REPRESENTATION ABOUT THE
 * SUITABILITY OF THESE LICENSED DELIVERABLES FOR ANY PURPOSE.  IT IS
 * PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF ANY KIND.
 * NVIDIA DISCLAIMS ALL WARRANTIES WITH REGARD TO THESE LICENSED
 * DELIVERABLES, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY,
 * NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * NOTWITHSTANDING ANY TERMS OR CONDITIONS TO THE CONTRARY IN THE
 * LICENSE AGREEMENT, IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY
 * SPECIAL, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OR ANY
 * DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
 * WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
 * ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
 * OF THESE LICENSED DELIVERABLES.
 *
 * U.S. Government End Users.  These Licensed Deliverables are a
 * "commercial item" as that term is defined at 48 C.F.R. 2.101 (OCT
 * 1995), consisting of "commercial computer software" and "commercial
 * computer software documentation" as such terms are used in 48
 * C.F.R. 12.212 (SEPT 1995) and is provided to the U.S. Government
 * only as a commercial end item.  Consistent with 48 C.F.R.12.212 and
 * 48 C.F.R. 227.7202-1 through 227.7202-4 (JUNE 1995), all
 * U.S. Government End Users acquire the Licensed Deliverables with
 * only those rights set forth herein.
 *
 * Any use of the Licensed Deliverables in individual and commercial
 * software must include, in the user documentation and internal
 * comments to the code, the above Disclaimer and U.S. Government End
 * Users Notice.
 */
#include "bicgstab/bicgstab.hpp"

#include <assert.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <stdio.h>   // fopen
#include <stdlib.h>  // EXIT_FAILURE
#include <string.h>  // strtok

/// This code allocates. The caller must free.
void make_test_matrix(int* n_out, int** row_offsets_out, int** columns_out, double** values_out) {
    int grid = 700;  // grid resolution

    int n = grid * grid;
    *n_out = n;
    // vertices have 5 neighbors,
    // but each vertex on the boundary loses 1. corners lose 2.
    int nnz = 5 * n - 4 * grid;

    printf(
        "Creating 5-point time-dependent advection-diffusion matrix.\n"
        " grid size: %d x %d\n"
        " matrix rows:   %d\n"
        " matrix cols:   %d\n"
        " nnz:         %d\n",
        grid, grid, n, n, nnz);

    int* row_offsets = *row_offsets_out = (int*)malloc((n + 1) * sizeof(int));
    int* columns = *columns_out = (int*)malloc(nnz * sizeof(int));
    double* values = *values_out = (double*)malloc(nnz * sizeof(double));
    assert(row_offsets);
    assert(columns);
    assert(values);

    double mass = 0.3;  // extra diagonal/mass term
    double ux = 0.3;    // advection velocity
    double uy = 0.2;

    int it = 0;  // next unused index into `columns`/`values`

#define INSERT(u, v, x)                                     \
    if (0 <= (u) && (u) < grid && 0 <= (v) && (v) < grid) { \
        columns[it] = ((u)*grid + (v));                     \
        values[it] = (x);                                   \
        ++it;                                               \
    }

    int row = 0;
    row_offsets[row] = 0;
    for (int i = 0; i < grid; ++i) {
        for (int j = 0; j < grid; ++j) {
            // Upwinding, so 'ux' and 'uy' only affect the -i and -j directions.
            INSERT(i - 1, j, -1.0 - ux);
            INSERT(i, j - 1, -1.0 - uy);
            INSERT(i, j, 4.0 + mass + ux + uy);
            INSERT(i, j + 1, -1.0);
            INSERT(i + 1, j, -1.0);
            row_offsets[++row] = it;
        }
    }
    assert(it == nnz);
#undef INSERT
}

//==============================================================================

int gpu_BiCGStab(cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle, int m,
                 cusparseSpMatDescr_t matA, cusparseSpMatDescr_t matM_lower,
                 cusparseSpMatDescr_t matM_upper, Vec d_B, Vec d_X, Vec d_R0, Vec d_R, Vec d_P,
                 Vec d_P_aux, Vec d_S, Vec d_S_aux, Vec d_V, Vec d_T, Vec d_tmp, void* d_bufferMV,
                 int maxIterations, double tolerance) {
    const double zero = 0.0;
    const double one = 1.0;
    const double minus_one = -1.0;
    //--------------------------------------------------------------------------
    // Create opaque data structures that holds analysis data between calls
    double coeff_tmp;
    size_t bufferSizeL, bufferSizeU;
    void *d_bufferL, *d_bufferU;
    cusparseSpSVDescr_t spsvDescrL, spsvDescrU;
    CHECK_CUSPARSE(cusparseSpSV_createDescr(&spsvDescrL))
    CHECK_CUSPARSE(cusparseSpSV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                           &coeff_tmp, matM_lower, d_P.vec, d_tmp.vec, CUDA_R_64F,
                                           CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, &bufferSizeL))
    CHECK_CUDA(cudaMalloc(&d_bufferL, bufferSizeL))
    CHECK_CUSPARSE(cusparseSpSV_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                         &coeff_tmp, matM_lower, d_P.vec, d_tmp.vec, CUDA_R_64F,
                                         CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL, d_bufferL))

    // Calculate UPPER buffersize
    CHECK_CUSPARSE(cusparseSpSV_createDescr(&spsvDescrU))
    CHECK_CUSPARSE(cusparseSpSV_bufferSize(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &coeff_tmp, matM_upper, d_tmp.vec,
        d_P_aux.vec, CUDA_R_64F, CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrU, &bufferSizeU))
    CHECK_CUDA(cudaMalloc(&d_bufferU, bufferSizeU))
    CHECK_CUSPARSE(cusparseSpSV_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                         &coeff_tmp, matM_upper, d_tmp.vec, d_P_aux.vec, CUDA_R_64F,
                                         CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrU, d_bufferU))
    //--------------------------------------------------------------------------
    // ### 1 ### R0 = b - A * X0 (using initial guess in X)
    //    (a) copy b in R0
    CHECK_CUDA(cudaMemcpy(d_R0.ptr, d_B.ptr, m * sizeof(double), cudaMemcpyDeviceToDevice))
    //    (b) compute R = -A * X0 + R
    CHECK_CUSPARSE(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &minus_one, matA,
                                d_X.vec, &one, d_R0.vec, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT,
                                d_bufferMV))
    //--------------------------------------------------------------------------
    double alpha, delta, delta_prev, omega;
    CHECK_CUBLAS(cublasDdot(cublasHandle, m, d_R0.ptr, 1, d_R0.ptr, 1, &delta))
    delta_prev = delta;
    // R = R0
    CHECK_CUDA(cudaMemcpy(d_R.ptr, d_R0.ptr, m * sizeof(double), cudaMemcpyDeviceToDevice))
    //--------------------------------------------------------------------------
    // nrm_R0 = ||R||
    double nrm_R;
    CHECK_CUBLAS(cublasDnrm2(cublasHandle, m, d_R0.ptr, 1, &nrm_R))
    double threshold = tolerance * nrm_R;
    printf("  Initial Residual: Norm %e' threshold %e\n", nrm_R, threshold);
    //--------------------------------------------------------------------------
    // ### 2 ### repeat until convergence based on max iterations and
    //           and relative residual
    for (int i = 1; i <= maxIterations; i++) {
        printf("  Iteration = %d; Error Norm = %e\n", i, nrm_R);
        //----------------------------------------------------------------------
        // ### 4, 7 ### P_i = R_i
        CHECK_CUDA(cudaMemcpy(d_P.ptr, d_R.ptr, m * sizeof(double), cudaMemcpyDeviceToDevice))
        if (i > 1) {
            //------------------------------------------------------------------
            // ### 6 ### beta = (delta_i / delta_i-1) * (alpha / omega_i-1)
            //    (a) delta_i = (R'_0, R_i-1)
            CHECK_CUBLAS(cublasDdot(cublasHandle, m, d_R0.ptr, 1, d_R.ptr, 1, &delta))
            //    (b) beta = (delta_i / delta_i-1) * (alpha / omega_i-1);
            double beta = (delta / delta_prev) * (alpha / omega);
            delta_prev = delta;
            //------------------------------------------------------------------
            // ### 7 ### P = R + beta * (P - omega * V)
            //    (a) P = - omega * V + P
            double minus_omega = -omega;
            CHECK_CUBLAS(cublasDaxpy(cublasHandle, m, &minus_omega, d_V.ptr, 1, d_P.ptr, 1))
            //    (b) P = beta * P
            CHECK_CUBLAS(cublasDscal(cublasHandle, m, &beta, d_P.ptr, 1))
            //    (c) P = R + P
            CHECK_CUBLAS(cublasDaxpy(cublasHandle, m, &one, d_R.ptr, 1, d_P.ptr, 1))
        }
        //----------------------------------------------------------------------
        // ### 9 ### P_aux = M_U^-1 M_L^-1 P_i
        //    (a) M_L^-1 P_i => tmp    (triangular solver)
        CHECK_CUDA(cudaMemset(d_tmp.ptr, 0x0, m * sizeof(double)))
        CHECK_CUDA(cudaMemset(d_P_aux.ptr, 0x0, m * sizeof(double)))
        CHECK_CUSPARSE(cusparseSpSV_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                          matM_lower, d_P.vec, d_tmp.vec, CUDA_R_64F,
                                          CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL))
        //    (b) M_U^-1 tmp => P_aux    (triangular solver)
        CHECK_CUSPARSE(cusparseSpSV_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                          matM_upper, d_tmp.vec, d_P_aux.vec, CUDA_R_64F,
                                          CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrU))
        //----------------------------------------------------------------------
        // ### 10 ### alpha = (R'0, R_i-1) / (R'0, A * P_aux)
        //    (a) V = A * P_aux
        CHECK_CUSPARSE(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, matA,
                                    d_P_aux.vec, &zero, d_V.vec, CUDA_R_64F,
                                    CUSPARSE_SPMV_ALG_DEFAULT, d_bufferMV))
        //    (b) denominator = R'0 * V
        double denominator;
        CHECK_CUBLAS(cublasDdot(cublasHandle, m, d_R0.ptr, 1, d_V.ptr, 1, &denominator))
        alpha = delta / denominator;
        PRINT_INFO(delta)
        PRINT_INFO(alpha)
        //----------------------------------------------------------------------
        // ### 11 ###  X_i = X_i-1 + alpha * P_aux
        CHECK_CUBLAS(cublasDaxpy(cublasHandle, m, &alpha, d_P_aux.ptr, 1, d_X.ptr, 1))
        //----------------------------------------------------------------------
        // ### 12 ###  S = R_i-1 - alpha * (A * P_aux)
        //    (a) S = R_i-1
        CHECK_CUDA(cudaMemcpy(d_S.ptr, d_R.ptr, m * sizeof(double), cudaMemcpyDeviceToDevice))
        //    (b) S = -alpha * V + R_i-1
        double minus_alpha = -alpha;
        CHECK_CUBLAS(cublasDaxpy(cublasHandle, m, &minus_alpha, d_V.ptr, 1, d_S.ptr, 1))
        //----------------------------------------------------------------------
        // ### 13 ###  check ||S|| < threshold
        double nrm_S;
        CHECK_CUBLAS(cublasDnrm2(cublasHandle, m, d_S.ptr, 1, &nrm_S))
        PRINT_INFO(nrm_S)
        if (nrm_S < threshold) break;
        //----------------------------------------------------------------------
        // ### 14 ### S_aux = M_U^-1 M_L^-1 S
        //    (a) M_L^-1 S => tmp    (triangular solver)
        cudaMemset(d_tmp.ptr, 0x0, m * sizeof(double));
        cudaMemset(d_S_aux.ptr, 0x0, m * sizeof(double));
        CHECK_CUSPARSE(cusparseSpSV_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                          matM_lower, d_S.vec, d_tmp.vec, CUDA_R_64F,
                                          CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrL))
        //    (b) M_U^-1 tmp => S_aux    (triangular solver)
        CHECK_CUSPARSE(cusparseSpSV_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                          matM_upper, d_tmp.vec, d_S_aux.vec, CUDA_R_64F,
                                          CUSPARSE_SPSV_ALG_DEFAULT, spsvDescrU))
        //----------------------------------------------------------------------
        // ### 15 ### omega = (A * S_aux, s) / (A * S_aux, A * S_aux)
        //    (a) T = A * S_aux
        CHECK_CUSPARSE(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, matA,
                                    d_S_aux.vec, &zero, d_T.vec, CUDA_R_64F,
                                    CUSPARSE_SPMV_ALG_DEFAULT, d_bufferMV))
        //    (b) omega_num = (A * S_aux, s)
        double omega_num, omega_den;
        CHECK_CUBLAS(cublasDdot(cublasHandle, m, d_T.ptr, 1, d_S.ptr, 1, &omega_num))
        //    (c) omega_den = (A * S_aux, A * S_aux)
        CHECK_CUBLAS(cublasDdot(cublasHandle, m, d_T.ptr, 1, d_T.ptr, 1, &omega_den))
        //    (d) omega = omega_num / omega_den
        omega = omega_num / omega_den;
        PRINT_INFO(omega)
        // ---------------------------------------------------------------------
        // ### 16 ### omega = X_i = X_i-1 + alpha * P_aux + omega * S_aux
        //    (a) X_i has been updated with h = X_i-1 + alpha * P_aux
        //        X_i = omega * S_aux + X_i
        CHECK_CUBLAS(cublasDaxpy(cublasHandle, m, &omega, d_S_aux.ptr, 1, d_X.ptr, 1))
        // ---------------------------------------------------------------------
        // ### 17 ###  R_i+1 = S - omega * (A * S_aux)
        //    (a) copy S in R
        CHECK_CUDA(cudaMemcpy(d_R.ptr, d_S.ptr, m * sizeof(double), cudaMemcpyDeviceToDevice))
        //    (a) R_i+1 = -omega * T + R
        double minus_omega = -omega;
        CHECK_CUBLAS(cublasDaxpy(cublasHandle, m, &minus_omega, d_T.ptr, 1, d_R.ptr, 1))
        // ---------------------------------------------------------------------
        // ### 18 ###  check ||R_i|| < threshold
        CHECK_CUBLAS(cublasDnrm2(cublasHandle, m, d_R.ptr, 1, &nrm_R))
        PRINT_INFO(nrm_R)
        if (nrm_R < threshold) break;
    }
    //--------------------------------------------------------------------------
    printf("Check Solution\n");  // ||R = b - A * X||
    //    (a) copy b in R
    CHECK_CUDA(cudaMemcpy(d_R.ptr, d_B.ptr, m * sizeof(double), cudaMemcpyDeviceToDevice))
    // R = -A * X + R
    CHECK_CUSPARSE(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &minus_one, matA,
                                d_X.vec, &one, d_R.vec, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT,
                                d_bufferMV))
    // check ||R||
    CHECK_CUBLAS(cublasDnrm2(cublasHandle, m, d_R.ptr, 1, &nrm_R))
    printf("Final error norm = %e\n", nrm_R);
    //--------------------------------------------------------------------------
    CHECK_CUSPARSE(cusparseSpSV_destroyDescr(spsvDescrL))
    CHECK_CUSPARSE(cusparseSpSV_destroyDescr(spsvDescrU))
    CHECK_CUDA(cudaFree(d_bufferL))
    CHECK_CUDA(cudaFree(d_bufferU))
    return EXIT_SUCCESS;
}