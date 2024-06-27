#include "../include/qmcp-solver/linear_programming_solver.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

#include "bam-api/read.hpp"
#include "logging/log.hpp"

void qmcp::LinearProgrammingSolver::make_matrix(int32_t* rows_out, int32_t** row_offsets_out,
                                                int32_t** columns_out, double** values_out) {
    bam_api::ReadIndex read_count = input_sequence_.get_reads_count();
    uint64_t rows = input_sequence_.ref_genome_length + read_count;
    *rows_out = rows;
    uint64_t nnz =
        read_count +  // we have I matrix of size readCount x readCount
        input_sequence_
            .ref_genome_length;  // and additional |G| columns with 1's starting from |R| index
    for (uint32_t i = 0; i < read_count; ++i) {
        nnz += input_sequence_.end_inds[i] - input_sequence_.start_inds[i] + 1;
    }

    int* row_offsets = *row_offsets_out = static_cast<int*>(malloc((rows + 1) * sizeof(int)));
    int* columns = *columns_out = (int*)malloc(nnz * sizeof(int));
    double* values = *values_out = (double*)malloc(nnz * sizeof(double));

    uint32_t value_ind = 0;

    // create Identity matrix
    for (uint32_t identity_row = 0; identity_row < read_count; ++identity_row) {
        row_offsets[identity_row] = value_ind;
        values[value_ind] = 1;
        columns[value_ind] = identity_row;
        value_ind++;
    }

    // create original matrix A
    for (uint32_t ref_ind_it = 0; ref_ind_it < input_sequence_.ref_genome_length; ++ref_ind_it) {
        row_offsets[read_count + ref_ind_it] = value_ind;
        for (uint32_t read_it = 0; read_it < input_sequence_.ids.size(); ++read_it) {
            if (input_sequence_.start_inds[read_it] <= ref_ind_it &&
                input_sequence_.end_inds[read_it] >= ref_ind_it) {
                values[value_ind] = -1;
                columns[value_ind] = read_it;
                // std::cout << "TERAZ KOLUMNA bedzie "
                value_ind++;
            }
        }
        // add |G| columns
        values[value_ind] = INT32_MAX;
        columns[value_ind] = ref_ind_it + read_count;
        value_ind++;
    }
    row_offsets[rows] = value_ind;
    LOG_WITH_LEVEL(logging::LogLevel::DEBUG) << "nnz: " << nnz << ", last offset: " << value_ind;
}

std::vector<double> qmcp::LinearProgrammingSolver::create_b_vector(uint32_t M) {
    bam_api::ReadIndex readCount = input_sequence_.get_reads_count();

    std::vector<double> b(input_sequence_.ref_genome_length + readCount, 0);

    for (uint32_t i = 0; i < readCount; i++) {
        b[i] = -1;
    }

    for (uint32_t i = 0; i < readCount; ++i) {
        for (uint32_t j = input_sequence_.start_inds[i]; j <= input_sequence_.end_inds[i]; ++j) {
            ++b[j + readCount];
        }
    }

    // cap nucleotides with more reads than M to M
    for (uint32_t i = 0; i < input_sequence_.ref_genome_length; ++i) {
        if (b[i + readCount] > M) b[i + readCount] = M;
    }

    return b;
}

std::vector<double> qmcp::LinearProgrammingSolver::process_bicgstab(int m, int* h_A_rows,
                                                                    int* h_A_columns,
                                                                    double* h_A_values,
                                                                    std::vector<double> h_B,
                                                                    std::vector<double> h_X) {
    // code from main.cpp
    const int maxIterations = 100;
    const double tolerance = 0.0000000001;

    int base = 0;

    int num_offsets = m + 1;
    int nnz = h_A_rows[m];
    //--------------------------------------------------------------------------
    // ### Device memory management ###
    int *d_A_rows, *d_A_columns;
    double *d_A_values, *d_M_values;
    Vec d_B, d_X, d_R, d_R0, d_P, d_P_aux, d_S, d_S_aux, d_V, d_T, d_tmp;

    // allocate device memory for CSR matrices
    CHECK_CUDA(cudaMalloc((void**)&d_A_rows, num_offsets * sizeof(int)))
    CHECK_CUDA(cudaMalloc((void**)&d_A_columns, nnz * sizeof(int)))
    CHECK_CUDA(cudaMalloc((void**)&d_A_values, nnz * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_M_values, nnz * sizeof(double)))

    CHECK_CUDA(cudaMalloc((void**)&d_B.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_X.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_R.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_R0.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_P.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_P_aux.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_S.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_S_aux.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_V.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_T.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_tmp.ptr, m * sizeof(double)))

    // copy the CSR matrices and vectors into device memory
    CHECK_CUDA(cudaMemcpy(d_A_rows, h_A_rows, num_offsets * sizeof(int), cudaMemcpyHostToDevice))
    CHECK_CUDA(cudaMemcpy(d_A_columns, h_A_columns, nnz * sizeof(int), cudaMemcpyHostToDevice))
    CHECK_CUDA(cudaMemcpy(d_A_values, h_A_values, nnz * sizeof(double), cudaMemcpyHostToDevice))
    CHECK_CUDA(cudaMemcpy(d_M_values, h_A_values, nnz * sizeof(double), cudaMemcpyHostToDevice))
    CHECK_CUDA(cudaMemcpy(d_X.ptr, h_X.data(), m * sizeof(double), cudaMemcpyHostToDevice))
    CHECK_CUDA(cudaMemcpy(d_B.ptr, h_B.data(), m * sizeof(double), cudaMemcpyHostToDevice))
    //--------------------------------------------------------------------------
    // ### cuSPARSE Handle and descriptors initialization ###
    // create the test matrix on the host
    cublasHandle_t cublasHandle = NULL;
    cusparseHandle_t cusparseHandle = NULL;
    CHECK_CUBLAS(cublasCreate(&cublasHandle))
    CHECK_CUSPARSE(cusparseCreate(&cusparseHandle))
    // Create dense vectors
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_B.vec, m, d_B.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_X.vec, m, d_X.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_R.vec, m, d_R.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_R0.vec, m, d_R0.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_P.vec, m, d_P.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_P_aux.vec, m, d_P_aux.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_S.vec, m, d_S.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_S_aux.vec, m, d_S_aux.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_V.vec, m, d_V.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_T.vec, m, d_T.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_tmp.vec, m, d_tmp.ptr, CUDA_R_64F))

    cusparseIndexBase_t baseIdx = CUSPARSE_INDEX_BASE_ZERO;
    // IMPORTANT: Upper/Lower triangular decompositions of A
    //            (matM_lower, matM_upper) must use two distinct descriptors
    cusparseSpMatDescr_t matA, matM_lower, matM_upper;
    cusparseMatDescr_t matLU;
    int* d_M_rows = d_A_rows;
    int* d_M_columns = d_A_columns;
    cusparseFillMode_t fill_lower = CUSPARSE_FILL_MODE_LOWER;
    cusparseDiagType_t diag_unit = CUSPARSE_DIAG_TYPE_UNIT;
    cusparseFillMode_t fill_upper = CUSPARSE_FILL_MODE_UPPER;
    cusparseDiagType_t diag_non_unit = CUSPARSE_DIAG_TYPE_NON_UNIT;
    // A
    CHECK_CUSPARSE(cusparseCreateCsr(&matA, m, m, nnz, d_A_rows, d_A_columns, d_A_values,
                                     CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, baseIdx, CUDA_R_64F))
    // M_lower
    CHECK_CUSPARSE(cusparseCreateCsr(&matM_lower, m, m, nnz, d_M_rows, d_M_columns, d_M_values,
                                     CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, baseIdx, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseSpMatSetAttribute(matM_lower, CUSPARSE_SPMAT_FILL_MODE, &fill_lower,
                                             sizeof(fill_lower)))
    CHECK_CUSPARSE(cusparseSpMatSetAttribute(matM_lower, CUSPARSE_SPMAT_DIAG_TYPE, &diag_unit,
                                             sizeof(diag_unit)))
    // M_upper
    CHECK_CUSPARSE(cusparseCreateCsr(&matM_upper, m, m, nnz, d_M_rows, d_M_columns, d_M_values,
                                     CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, baseIdx, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseSpMatSetAttribute(matM_upper, CUSPARSE_SPMAT_FILL_MODE, &fill_upper,
                                             sizeof(fill_upper)))
    CHECK_CUSPARSE(cusparseSpMatSetAttribute(matM_upper, CUSPARSE_SPMAT_DIAG_TYPE, &diag_non_unit,
                                             sizeof(diag_non_unit)))
    //--------------------------------------------------------------------------
    // ### Preparation ### b = A * X
    const double alpha = 0.75;
    size_t bufferSizeMV;
    void* d_bufferMV;
    double beta = 0.0;
    CHECK_CUSPARSE(cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha,
                                           matA, d_X.vec, &beta, d_B.vec, CUDA_R_64F,
                                           CUSPARSE_SPMV_ALG_DEFAULT, &bufferSizeMV))
    CHECK_CUDA(cudaMalloc(&d_bufferMV, bufferSizeMV))

    CHECK_CUSPARSE(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA,
                                d_X.vec, &beta, d_B.vec, CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT,
                                d_bufferMV))
    // X0 = 0
    CHECK_CUDA(cudaMemset(d_X.ptr, 0x0, m * sizeof(double)))
    //--------------------------------------------------------------------------
    // Perform Incomplete-LU factorization of A (csrilu0) -> M_lower, M_upper
    csrilu02Info_t infoM = NULL;
    int bufferSizeLU = 0;
    void* d_bufferLU;
    CHECK_CUSPARSE(cusparseCreateMatDescr(&matLU))
    CHECK_CUSPARSE(cusparseSetMatType(matLU, CUSPARSE_MATRIX_TYPE_GENERAL))
    CHECK_CUSPARSE(cusparseSetMatIndexBase(matLU, baseIdx))
    CHECK_CUSPARSE(cusparseCreateCsrilu02Info(&infoM))

    CHECK_CUSPARSE(cusparseDcsrilu02_bufferSize(cusparseHandle, m, nnz, matLU, d_M_values, d_A_rows,
                                                d_A_columns, infoM, &bufferSizeLU))
    CHECK_CUDA(cudaMalloc(&d_bufferLU, bufferSizeLU))
    CHECK_CUSPARSE(cusparseDcsrilu02_analysis(cusparseHandle, m, nnz, matLU, d_M_values, d_A_rows,
                                              d_A_columns, infoM, CUSPARSE_SOLVE_POLICY_USE_LEVEL,
                                              d_bufferLU))
    int structural_zero = 0;
    CHECK_CUSPARSE(cusparseXcsrilu02_zeroPivot(cusparseHandle, infoM, &structural_zero))

    //  M = L * U
    CHECK_CUSPARSE(cusparseDcsrilu02(cusparseHandle, m, nnz, matLU, d_M_values, d_A_rows,
                                     d_A_columns, infoM, CUSPARSE_SOLVE_POLICY_USE_LEVEL,
                                     d_bufferLU))
    // Find numerical zero
    int numerical_zero;
    CHECK_CUSPARSE(cusparseXcsrilu02_zeroPivot(cusparseHandle, infoM, &numerical_zero))

    CHECK_CUSPARSE(cusparseDestroyCsrilu02Info(infoM))
    CHECK_CUSPARSE(cusparseDestroyMatDescr(matLU))
    CHECK_CUDA(cudaFree(d_bufferLU))
    //--------------------------------------------------------------------------
    // ### Run BiCGStab computation ###
    // printf("BiCGStab loop:\n");
    gpu_BiCGStab(cublasHandle, cusparseHandle, m, matA, matM_lower, matM_upper, d_B, d_X, d_R0, d_R,
                 d_P, d_P_aux, d_S, d_S_aux, d_V, d_T, d_tmp, d_bufferMV, maxIterations, tolerance);
    //--------------------------------------------------------------------------
    // ### Free resources ###
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_B.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_X.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_R.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_R0.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_P.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_P_aux.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_S.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_S_aux.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_V.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_T.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_tmp.vec))
    CHECK_CUSPARSE(cusparseDestroySpMat(matA))
    CHECK_CUSPARSE(cusparseDestroySpMat(matM_lower))
    CHECK_CUSPARSE(cusparseDestroySpMat(matM_upper))
    CHECK_CUSPARSE(cusparseDestroy(cusparseHandle))
    CHECK_CUBLAS(cublasDestroy(cublasHandle))

    free(h_A_rows);
    free(h_A_columns);
    free(h_A_values);

    CHECK_CUDA(cudaFree(d_X.ptr))
    CHECK_CUDA(cudaFree(d_B.ptr))
    CHECK_CUDA(cudaFree(d_R.ptr))
    CHECK_CUDA(cudaFree(d_R0.ptr))
    CHECK_CUDA(cudaFree(d_P.ptr))
    CHECK_CUDA(cudaFree(d_P_aux.ptr))
    CHECK_CUDA(cudaFree(d_S.ptr))
    CHECK_CUDA(cudaFree(d_S_aux.ptr))
    CHECK_CUDA(cudaFree(d_V.ptr))
    CHECK_CUDA(cudaFree(d_T.ptr))
    CHECK_CUDA(cudaFree(d_tmp.ptr))
    CHECK_CUDA(cudaFree(d_A_values))
    CHECK_CUDA(cudaFree(d_A_columns))
    CHECK_CUDA(cudaFree(d_A_rows))
    CHECK_CUDA(cudaFree(d_M_values))
    CHECK_CUDA(cudaFree(d_bufferMV))
    return h_X;
}

std::unique_ptr<qmcp::Solution> qmcp::LinearProgrammingSolver::solve(uint32_t max_coverage,
                                                                     bam_api::BamApi& bam_api) {
    input_sequence_ = bam_api.get_paired_reads_soa();

    int m = -1;
    int* h_A_rows = NULL;
    int* h_A_columns = NULL;
    double* h_A_values = NULL;

    make_matrix(&m, &h_A_rows, &h_A_columns, &h_A_values);
    auto b = create_b_vector(max_coverage);

    auto x = process_bicgstab(m, h_A_rows, h_A_columns, h_A_values, b, std::vector<double>(m, 1));

    // return solution
    auto reduced_reads = std::make_unique<Solution>();

    for (bam_api::ReadIndex i = 0; i < input_sequence_.get_reads_count(); ++i) {
        if (x[i] > 0) {
            // LOG_WITH_LEVEL(logging::DEBUG) << "h_X[" << i << "]: " << h_X[i];
            reduced_reads->push_back(input_sequence_.ids[i]);
        }
    }
    return reduced_reads;
}
