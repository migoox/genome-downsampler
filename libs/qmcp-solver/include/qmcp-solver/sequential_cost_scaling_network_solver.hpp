#pragma once

#include <ortools/graph/min_cost_flow.h>

#include <filesystem>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"
#include "solver.hpp"

namespace qmcp {
class SequentialCostScalingNetworkSolver : public Solver {
   public:
    SequentialCostScalingNetworkSolver();
    explicit SequentialCostScalingNetworkSolver(const std::filesystem::path& filepath,
                                                uint32_t min_seq_length, uint32_t min_seq_mapq);

    void import_reads(const std::filesystem::path& filepath, uint32_t min_seq_length,
                      uint32_t min_seq_mapq) override;
    void solve(uint32_t max_coverage) override;
    void export_reads(const std::filesystem::path& filepath) override;

    std::vector<bam_api::ReadIndex> output_sequence();

   private:
    static void create_network_flow_graph(operations_research::SimpleMinCostFlow& min_cost_flow,
                                          const bam_api::AOSPairedReads& sequence, unsigned int M);
    static std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence,
                                                   unsigned int M);
    static std::vector<bam_api::ReadIndex> obtain_sequence(
        const bam_api::AOSPairedReads& sequence,
        const operations_research::SimpleMinCostFlow& min_cost_flow);

    static void add_pairs(std::vector<bam_api::ReadIndex>& reduced_reads,
                          const std::vector<bool>& mapped_reads,
                          const std::vector<std::optional<bam_api::ReadIndex>>& read_pair_map);

    std::filesystem::path input_filepath_;
    bam_api::AOSPairedReads input_sequence_;
    std::vector<bam_api::ReadIndex> output_sequence_;
    bool is_data_loaded_ = false;
};
}  // namespace qmcp
