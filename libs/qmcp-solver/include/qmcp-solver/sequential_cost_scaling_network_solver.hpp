#ifndef QMCP_SEQUENTIAL_COST_SCALING_NETWORK_SOLVER_HPP
#define QMCP_SEQUENTIAL_COST_SCALING_NETWORK_SOLVER_HPP
#include <filesystem>
#include <iostream>
#include "solver.hpp"
#include "../../bam-api/include/bam-api/bam_paired_reads.hpp"
#include "../../bam-api/include/bam-api/bam_api.hpp"
#include <ortools/graph/min_cost_flow.h>


namespace qmcp {
    
class SequentialCostScalingNetworkSolver : public Solver {
   public:
    SequentialCostScalingNetworkSolver(unsigned int M,const std::filesystem::path& filepath) : M_(M) {
        input_sequence_ = bam_api::BamApi::read_bam_aos(filepath);
        std::cout<<"read bam\n";
    }
    void solve() override;
    private:
    bam_api::AOSPairedReads input_sequence_;
    bam_api::AOSPairedReads output_sequence_;
    unsigned int M_;

    static void create_network_flow_graph(operations_research::SimpleMinCostFlow& min_cost_flow, const bam_api::AOSPairedReads& sequence, unsigned int M);
    
};

}  // namespace qmcp

#endif
