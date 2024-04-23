#ifndef QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#define QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#include <cstdio>
#include <filesystem>
#include "solver.hpp"
#include "../../bam-api/include/bam-api/bam_paired_reads.hpp"
#include "../../bam-api/include/bam-api/bam_api.hpp"



namespace qmcp {

class SequenceNetworkSolver : public Solver {
   public:
    SequenceNetworkSolver(unsigned int M,const std::filesystem::path& filepath) : M_(M) {
        printf("START\n");
        sequence_ = bam_api::BamApi::read_bam_aos(filepath);
        printf("KONIEC\n");
    }
    void solve() override;
    private:
    bam_api::AOSPairedReads sequence_;
    unsigned int M_;
    
};

}  // namespace qmcp

#endif
