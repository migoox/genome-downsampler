
#include "qmcp-solver/conjugate_gradient_solver.hpp"
#include <driver_types.h>

#include <algorithm>
#include <cstdint>
#include <iostream>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"
#include "cg/cg.hpp"
#include "logging/log.hpp"

void qmcp::ConjugateGradientSolver::import_reads(const std::filesystem::path& input,
                                                 uint32_t min_seq_length, uint32_t min_seq_mapq) {
    input_ = input;
    LOG_WITH_LEVEL(logging::kDebug) << "Import, min_len: " << min_seq_length
                                    << ", min_mapq: " << min_seq_mapq << ", input: " << input;

    paired_reads_ = bam_api::BamApi::read_bam_soa(input, min_seq_length, min_seq_mapq);

    LOG_WITH_LEVEL(logging::kInfo) << paired_reads_.ids.size() << " sequences has been imported!";
}

void qmcp::ConjugateGradientSolver::make_matrix(int32_t* n_out, int32_t** row_offsets_out,
                                                int32_t** columns_out, double** values_out) {
    uint64_t n = paired_reads_.ref_genome_length;
    *n_out = n;
    uint64_t nnz = 0;
    for (uint32_t i = 0; i < paired_reads_.ids.size(); ++i) {
        nnz += paired_reads_.end_inds[i] - paired_reads_.start_inds[i] + 1;
    }

    int* row_offsets = *row_offsets_out = (int*)malloc((n + 1) * sizeof(int));
    int* columns = *columns_out = (int*)malloc(nnz * sizeof(int));
    double* values = *values_out = (double*)malloc(nnz * sizeof(double));

    uint32_t value_ind = 0;
    for(uint32_t ref_ind_it = 0; ref_ind_it < paired_reads_.ref_genome_length; ++ref_ind_it) {
      row_offsets[ref_ind_it] = value_ind;
      for(uint32_t read_it = 0; read_it < paired_reads_.ids.size(); ++read_it) {
          if(paired_reads_.start_inds[read_it] <= ref_ind_it && paired_reads_.end_inds[read_it] >= ref_ind_it) {
            values[value_ind] = 1;
            columns[value_ind] = read_it;
            value_ind++;
          }
      }
    }

    row_offsets[n] = value_ind;

    LOG_WITH_LEVEL(logging::kDebug) << "nnz: " << nnz << ", last offset: " << value_ind;
}

std::vector<double> qmcp::ConjugateGradientSolver::create_b_vector(uint32_t M) {
    std::vector<double> x(paired_reads_.ref_genome_length, 0);

    for (uint32_t i = 0; i < paired_reads_.ids.size(); ++i) {
        for (uint32_t j = paired_reads_.start_inds[i]; j <= paired_reads_.end_inds[i]; ++j) {
            ++x[j];
        }
    }

    // cap nucleotides with more reads than M to M
    for (uint32_t i = 0; i < paired_reads_.ref_genome_length; ++i) {
        if (x[i] > M) x[i] = M;
    }

    return x;
}

void qmcp::ConjugateGradientSolver::solve(uint32_t max_coverage) {
    LOG_WITH_LEVEL(logging::kDebug) << "Solve (max_coverage set to " << max_coverage << ")";

    const int32_t maxIterations = 10000;
    const double tolerance = 1e-8f;

    int32_t base = 0;
    int32_t m = -1;
    int32_t* h_A_rows = NULL;
    int32_t* h_A_columns = NULL;
    double* h_A_values = NULL;
    make_matrix(&m, &h_A_rows, &h_A_columns, &h_A_values);
    int32_t num_offsets = m + 1;
    int32_t nnz = h_A_rows[m];

    std::vector<double> h_B = create_b_vector(max_coverage);
    //--------------------------------------------------------------------------
    // ### Device memory management ###
    int32_t *d_A_rows, *d_A_columns;
    double *d_A_values, *d_L_values;
    Vec d_B, d_X, d_R, d_R_aux, d_P, d_T, d_tmp;

    // allocate device memory for CSR matrices
    CHECK_CUDA(cudaMalloc((void**)&d_A_rows, num_offsets * sizeof(int32_t)))
    CHECK_CUDA(cudaMalloc((void**)&d_A_columns, nnz * sizeof(int32_t)))
    CHECK_CUDA(cudaMalloc((void**)&d_A_values, nnz * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_L_values, nnz * sizeof(double)))

    CHECK_CUDA(cudaMalloc((void**)&d_B.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_X.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_R.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_R_aux.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_P.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_T.ptr, m * sizeof(double)))
    CHECK_CUDA(cudaMalloc((void**)&d_tmp.ptr, m * sizeof(double)))

    // copy the CSR matrices and vectors into device memory
    CHECK_CUDA(cudaMemcpy(d_A_rows, h_A_rows, num_offsets * sizeof(int32_t), cudaMemcpyHostToDevice))
    CHECK_CUDA(cudaMemcpy(d_A_columns, h_A_columns, nnz * sizeof(int32_t), cudaMemcpyHostToDevice))
    CHECK_CUDA(cudaMemcpy(d_A_values, h_A_values, nnz * sizeof(double), cudaMemcpyHostToDevice))
    CHECK_CUDA(cudaMemcpy(d_L_values, h_A_values, nnz * sizeof(double), cudaMemcpyHostToDevice))
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
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_R_aux.vec, m, d_R_aux.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_P.vec, m, d_P.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_T.vec, m, d_T.ptr, CUDA_R_64F))
    CHECK_CUSPARSE(cusparseCreateDnVec(&d_tmp.vec, m, d_tmp.ptr, CUDA_R_64F))

    cusparseIndexBase_t baseIdx = CUSPARSE_INDEX_BASE_ZERO;
    cusparseSpMatDescr_t matA, matL;
    int32_t* d_L_rows = d_A_rows;
    int32_t* d_L_columns = d_A_columns;
    cusparseFillMode_t fill_lower = CUSPARSE_FILL_MODE_LOWER;
    cusparseDiagType_t diag_non_unit = CUSPARSE_DIAG_TYPE_NON_UNIT;
    // A
    CHECK_CUSPARSE(cusparseCreateCsr(&matA, m, m, nnz, d_A_rows, d_A_columns, d_A_values,
                                     CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, baseIdx, CUDA_R_64F))
    // L
    CHECK_CUSPARSE(cusparseCreateCsr(&matL, m, m, nnz, d_L_rows, d_L_columns, d_L_values,
                                     CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, baseIdx, CUDA_R_64F))
    CHECK_CUSPARSE(
        cusparseSpMatSetAttribute(matL, CUSPARSE_SPMAT_FILL_MODE, &fill_lower, sizeof(fill_lower)))
    CHECK_CUSPARSE(cusparseSpMatSetAttribute(matL, CUSPARSE_SPMAT_DIAG_TYPE, &diag_non_unit,
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
    // Perform Incomplete-Cholesky factorization of A (csric0) -> L, L^T
    cusparseMatDescr_t descrM;
    csric02Info_t infoM = NULL;
    int32_t bufferSizeIC = 0;
    void* d_bufferIC;
    CHECK_CUSPARSE(cusparseCreateMatDescr(&descrM))
    CHECK_CUSPARSE(cusparseSetMatIndexBase(descrM, baseIdx))
    CHECK_CUSPARSE(cusparseSetMatType(descrM, CUSPARSE_MATRIX_TYPE_GENERAL))
    CHECK_CUSPARSE(cusparseSetMatFillMode(descrM, CUSPARSE_FILL_MODE_LOWER))
    CHECK_CUSPARSE(cusparseSetMatDiagType(descrM, CUSPARSE_DIAG_TYPE_NON_UNIT))
    CHECK_CUSPARSE(cusparseCreateCsric02Info(&infoM))

    CHECK_CUSPARSE(cusparseDcsric02_bufferSize(cusparseHandle, m, nnz, descrM, d_L_values, d_A_rows,
                                               d_A_columns, infoM, &bufferSizeIC))
    CHECK_CUDA(cudaMalloc(&d_bufferIC, bufferSizeIC))
    CHECK_CUSPARSE(cusparseDcsric02_analysis(cusparseHandle, m, nnz, descrM, d_L_values, d_A_rows,
                                             d_A_columns, infoM, CUSPARSE_SOLVE_POLICY_NO_LEVEL,
                                             d_bufferIC))
    int32_t structural_zero;
    CHECK_CUSPARSE(cusparseXcsric02_zeroPivot(cusparseHandle, infoM, &structural_zero))
    // M = L * L^T
    CHECK_CUSPARSE(cusparseDcsric02(cusparseHandle, m, nnz, descrM, d_L_values, d_A_rows,
                                    d_A_columns, infoM, CUSPARSE_SOLVE_POLICY_NO_LEVEL, d_bufferIC))
    // Find numerical zero
    int32_t numerical_zero;
    CHECK_CUSPARSE(cusparseXcsric02_zeroPivot(cusparseHandle, infoM, &numerical_zero))

    CHECK_CUSPARSE(cusparseDestroyCsric02Info(infoM))
    CHECK_CUSPARSE(cusparseDestroyMatDescr(descrM))
    CHECK_CUDA(cudaFree(d_bufferIC))
    //--------------------------------------------------------------------------
    // ### Run CG computation ###
    printf("CG loop:\n");
    gpu_CG(cublasHandle, cusparseHandle, m, matA, matL, d_B, d_X, d_R, d_R_aux, d_P, d_T, d_tmp,
           d_bufferMV, maxIterations, tolerance);


    std::vector<double> h_X;
    CHECK_CUDA(cudaMemcpy(h_X.data(), d_X.ptr, h_B.size() * sizeof(double), cudaMemcpyDeviceToHost))
    //--------------------------------------------------------------------------
    // ### Free resources ###
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_B.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_X.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_R.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_R_aux.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_P.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_T.vec))
    CHECK_CUSPARSE(cusparseDestroyDnVec(d_tmp.vec))
    CHECK_CUSPARSE(cusparseDestroySpMat(matA))
    CHECK_CUSPARSE(cusparseDestroySpMat(matL))
    CHECK_CUSPARSE(cusparseDestroy(cusparseHandle))
    CHECK_CUBLAS(cublasDestroy(cublasHandle))

    free(h_A_rows);
    free(h_A_columns);
    free(h_A_values);

    CHECK_CUDA(cudaFree(d_X.ptr))
    CHECK_CUDA(cudaFree(d_B.ptr))
    CHECK_CUDA(cudaFree(d_R.ptr))
    CHECK_CUDA(cudaFree(d_R_aux.ptr))
    CHECK_CUDA(cudaFree(d_P.ptr))
    CHECK_CUDA(cudaFree(d_T.ptr))
    CHECK_CUDA(cudaFree(d_tmp.ptr))
    CHECK_CUDA(cudaFree(d_A_values))
    CHECK_CUDA(cudaFree(d_A_columns))
    CHECK_CUDA(cudaFree(d_A_rows))
    CHECK_CUDA(cudaFree(d_L_values))
    CHECK_CUDA(cudaFree(d_bufferMV))

    for(bam_api::ReadIndex i = 0; i < paired_reads_.get_reads_count(); ++i) {
        if(h_X[i] > 0) {
            LOG_WITH_LEVEL(logging::kDebug) << "h_X[" << i << "]: " << h_X[i];
            solution_.push_back(paired_reads_.ids[i]);
        }
    }

    LOG_WITH_LEVEL(logging::kInfo) << "Solution have " << solution_.size() << " sequences!";
}


void qmcp::ConjugateGradientSolver::find_pairs(bool flag) {
  find_pairs_ = flag;
}

void qmcp::ConjugateGradientSolver::set_reads(const bam_api::SOAPairedReads& input_sequence) {
  paired_reads_ = input_sequence;
}

const std::vector<bam_api::BAMReadId>& qmcp::ConjugateGradientSolver::get_output() {
    return solution_;
}

std::vector<bam_api::BAMReadId> qmcp::ConjugateGradientSolver::output_sequence() {
    return solution_;
}

void qmcp::ConjugateGradientSolver::export_reads(const std::filesystem::path& output) {
    LOG_WITH_LEVEL(logging::kDebug) << "Output: " << output;

    uint32_t reads_written = 0;
    reads_written = bam_api::BamApi::write_bam(input_, output, solution_);

    LOG_WITH_LEVEL(logging::kInfo)
        << reads_written << " sequences has been written to file " << output;
}
