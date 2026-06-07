#include "../include/qmcp-solver/qmcp_cpu_cost_scaling_solver.hpp"

#include <memory>
#include <vector>

#include "bam-api/read.hpp"
#include "bam-api/region_api.hpp"
#include "logging/log.hpp"
#include "qmcp-solver/ortools_flow_types.hpp"
#include "qmcp-solver/solver.hpp"

std::unique_ptr<qmcp::Solution> qmcp::QmcpCpuCostScalingSolver::solve(uint32_t max_coverage,
                                                                      bam_api::RegionApi& reg_api) {
    input_sequence_ = reg_api.get_paired_reads_aos();

    operations_research::SimpleMinCostFlow min_cost_flow;

    create_network_flow_graph(min_cost_flow, input_sequence_, max_coverage);

    const int status = min_cost_flow.Solve();
    LOG_WITH_LEVEL(logging::DEBUG) << "flow == " << min_cost_flow.MaximumFlow()
                                   << ", optimal cost == " << min_cost_flow.OptimalCost();

    if (status != operations_research::MinCostFlow::OPTIMAL) {
        LOG_WITH_LEVEL(logging::ERROR)
            << "Solving the min cost flow problem failed. Solver status: " << status;
    }

    return obtain_sequence(input_sequence_, min_cost_flow);
}

void qmcp::QmcpCpuCostScalingSolver::create_network_flow_graph(
    operations_research::SimpleMinCostFlow& min_cost_flow, const bam_api::AOSPairedReads& sequence,
    uint32_t M) {
    const int capacity_upper_bound_multiplier = 100;
    const OrNode ref_len = to_node(sequence.ref_genome_length);
    const OrFlow backward_capacity =
        to_flow(static_cast<int64_t>(capacity_upper_bound_multiplier) * M);

    bam_api::ReadQuality max_quality = 0;
    for (const bam_api::Read& read : sequence.reads) {
        if (read.quality > max_quality) max_quality = read.quality;
    }

    for (const bam_api::Read& read : sequence.reads) {
        const OrCost cost = to_cost(max_quality - read.quality + 1);
        min_cost_flow.AddArcWithCapacityAndUnitCost(to_node(read.start_ind),
                                                    to_node(read.end_ind + 1), to_flow(1), cost);
    }

    for (OrNode i = 0; i < ref_len; ++i) {
        min_cost_flow.AddArcWithCapacityAndUnitCost(i + 1, i, backward_capacity, to_cost(0));
    }

    const std::vector<int> demand = create_demand_function(sequence, M);
    for (OrNode i = 0; i <= ref_len; ++i) {
        // ATTENTION!
        //  OR-Tools defines supply and demand inversely to our definition!
        //  So below it sets NodeSupply to -demand[i] where
        //  demand is function d from our understanding of a problem
        min_cost_flow.SetNodeSupply(i, to_signed_flow(-demand[static_cast<std::size_t>(i)]));
    }
}

std::vector<int> qmcp::QmcpCpuCostScalingSolver::create_b_function(
    const bam_api::AOSPairedReads& sequence, uint32_t M) {
    std::vector<int> b(sequence.ref_genome_length + 1, 0);

    for (const bam_api::Read& read : sequence.reads) {
        for (uint64_t j = read.start_ind; j <= read.end_ind; ++j) {
            ++b[j + 1];
        }
    }

    for (uint64_t i = 0; i <= sequence.ref_genome_length; ++i) {
        if (b[i] > static_cast<int>(M)) b[i] = static_cast<int>(M);
    }

    return b;
}

std::vector<int> qmcp::QmcpCpuCostScalingSolver::create_demand_function(
    const bam_api::AOSPairedReads& sequence, uint32_t M) {
    std::vector<int> b = create_b_function(sequence, M);

    const int b_1 = b[1];

    for (std::size_t i = 1; i < sequence.ref_genome_length; ++i) {
        b[i] -= b[i + 1];
    }

    b[0] = -b_1;

    return b;
}
std::unique_ptr<qmcp::Solution> qmcp::QmcpCpuCostScalingSolver::obtain_sequence(
    const bam_api::AOSPairedReads& sequence,
    const operations_research::SimpleMinCostFlow& min_cost_flow) {
    auto reduced_reads = std::make_unique<Solution>();

    for (bam_api::ReadIndex read_id = 0; read_id < sequence.get_reads_count(); ++read_id) {
        if (min_cost_flow.Flow(to_arc(read_id)) > 0) {
            reduced_reads->push_back(read_id);
        }
    }

    return reduced_reads;
}
