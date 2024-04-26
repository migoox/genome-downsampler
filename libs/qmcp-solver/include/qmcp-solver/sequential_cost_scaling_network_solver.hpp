#ifndef QMCP_SEQUENTIAL_COST_SCALING_NETWORK_SOLVER_HPP
#define QMCP_SEQUENTIAL_COST_SCALING_NETWORK_SOLVER_HPP
#include <cstdio>
#include <filesystem>
#include <iostream>
#include "solver.hpp"
#include "../../bam-api/include/bam-api/bam_paired_reads.hpp"
#include "../../bam-api/include/bam-api/bam_api.hpp"



namespace qmcp {

class SequentialCostScalingNetworkSolver : public Solver {
   public:
    SequentialCostScalingNetworkSolver(unsigned int M, const std::filesystem::path& filepath) : M_(M) {
        sequence_ = bam_api::BamApi::read_bam_aos(filepath);
        std::cout<<"read bam\n";
    }
    void solve() override;
    private:
    bam_api::AOSPairedReads sequence_;
    unsigned int M_;
    
};

}  // namespace qmcp

#endif
