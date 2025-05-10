#include "../include/qmcp-solver/qmcp_cpu_cost_scaling_solver.hpp"

#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/read.hpp"
#include "logging/log.hpp"
#include "qmcp-solver/solver.hpp"

std::unique_ptr<qmcp::Solution> qmcp::QmcpCpuCostScalingSolver::solve(uint32_t max_coverage,
                                                                      bam_api::BamApi& bam_api) {
    input_sequence_ = bam_api.get_paired_reads_aos();

    operations_research::SimpleMinCostFlow min_cost_flow;

    create_network_flow_graph(min_cost_flow, input_sequence_, max_coverage);

    int status = min_cost_flow.Solve();
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

    // 1. Find max quality
    bam_api::ReadQuality max_quality = 0;
    for (const bam_api::Read& read : sequence.reads) {
        if (read.quality > max_quality) max_quality = read.quality;
    }

    // 2. Create normal edges with calculated cost
    for (const bam_api::Read& read : sequence.reads) {
        int cost = static_cast<int>(max_quality - read.quality + 1);
        min_cost_flow.AddArcWithCapacityAndUnitCost(static_cast<int>(read.start_ind),
                                                    static_cast<int>(read.end_ind + 1), 1, cost);
    }

    // 3. Create backwards edges with cost 0
    //  it is safe to say that algorithm won't use more than c_u_b_m * M
    //  flow in one edge.
    for (int i = 0; i < sequence.ref_genome_length; ++i) {
        min_cost_flow.AddArcWithCapacityAndUnitCost(i + 1, i, capacity_upper_bound_multiplier * M,
                                                    0);
    }

    // 4. Add supply and demand (negative supply = demand)
    std::vector<int> demand = create_demand_function(sequence, M);
    for (int i = 0; i <= sequence.ref_genome_length; ++i) {
        // ATTENTION!
        //  OR-Tools defines supply and demand inversely to our definition!
        //  So below it sets NodeSupply to -demand[i] where
        //  demand is function d from our understanding of a problem
        min_cost_flow.SetNodeSupply(i, -demand[i]);
    }
}

std::vector<int> qmcp::QmcpCpuCostScalingSolver::create_b_function(
    const bam_api::AOSPairedReads& sequence, uint32_t M) {
    std::vector<int> b(sequence.ref_genome_length + 1, 0);

    for (const bam_api::Read& read : sequence.reads) {
        for (uint32_t j = read.start_ind; j <= read.end_ind; ++j) {
            ++b[j + 1];
        }
    }

    // cap nucleotides with more reads than M to M
    for (uint32_t i = 0; i <= sequence.ref_genome_length; ++i) {
        if (b[i] > M) b[i] = static_cast<int>(M);
    }

    return b;
}

std::vector<int> qmcp::QmcpCpuCostScalingSolver::create_demand_function(
    const bam_api::AOSPairedReads& sequence, uint32_t M) {
    std::vector<int> b = create_b_function(sequence, M);

    int b_1 = b[1];

    for (int i = 1; i < sequence.ref_genome_length; ++i) {
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
        if (min_cost_flow.Flow(static_cast<int>(read_id)) > 0) {
            reduced_reads->push_back(read_id);
        }
    }

    return reduced_reads;
}
