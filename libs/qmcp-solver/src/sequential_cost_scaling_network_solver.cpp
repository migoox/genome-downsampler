#include "../include/qmcp-solver/sequential_cost_scaling_network_solver.hpp"

#include <cstdint>
#include <iostream>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"

std::vector<bam_api::BAMReadId> qmcp::SequentialCostScalingNetworkSolver::solve(uint32_t max_coverage) {
    input_sequence_ = bam_api_.get_paired_reads_aos();

    operations_research::SimpleMinCostFlow min_cost_flow;

    create_network_flow_graph(min_cost_flow, input_sequence_, max_coverage);

    int status = min_cost_flow.Solve();

    if (status != operations_research::MinCostFlow::OPTIMAL) {
        LOG(INFO) << "Solving the min cost flow problem failed. Solver status: " << status;
    }
    
    return obtain_sequence(input_sequence_, min_cost_flow);
}

void qmcp::SequentialCostScalingNetworkSolver::create_network_flow_graph(
    operations_research::SimpleMinCostFlow& min_cost_flow, const bam_api::AOSPairedReads& sequence,
    unsigned int M) {
    const int capacity_upper_bound_multiplier = 100;
    // it is safe to say that algorithm wont use more than c_u_b_m * M
    // flow in one edge.

    // create normal edges
    for (int i = 0; i < sequence.reads.size(); ++i) {
        int weight =
            (sequence.reads[i].start_ind - sequence.reads[i].end_ind) * sequence.reads[i].quality;
        min_cost_flow.AddArcWithCapacityAndUnitCost(sequence.reads[i].start_ind - 1,
                                                    sequence.reads[i].end_ind, 1, weight);
    }

    // create backwards edge
    for (int i = 0; i < sequence.ref_genome_length; ++i) {
        min_cost_flow.AddArcWithCapacityAndUnitCost(i + 1, i, M * capacity_upper_bound_multiplier,
                                                    0);
    }

    min_cost_flow.AddArcWithCapacityAndUnitCost(0, sequence.ref_genome_length,
                                                M * capacity_upper_bound_multiplier, 0);

    // add supply and demand (negative supply = demand)
    std::vector<int> demand = create_demand_function(sequence, M);
    for (int i = 0; i < sequence.ref_genome_length; ++i) {
        min_cost_flow.SetNodeSupply(i, demand[i]);
    }
}

std::vector<int> create_b_function(const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int> b(sequence.ref_genome_length + 1, 0);

    for (unsigned int i = 0; i < sequence.reads.size(); ++i) {
        for (unsigned int j = sequence.reads[i].start_ind; j < sequence.reads[i].end_ind; ++j) {
            ++b[j];
        }
    }

    // cap nucleotides with more reads than M to M
    for (unsigned int i = 0; i < sequence.ref_genome_length + 1; ++i) {
        if (b[i] > M) b[i] = M;
    }
    return b;
}

std::vector<int> qmcp::SequentialCostScalingNetworkSolver::create_demand_function(
    const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int> b = create_b_function(sequence, M);

    int b_0 = b[0];

    for (int i = 0; i < sequence.ref_genome_length - 1; ++i) {
        b[i] = b[i] - b[i + 1];
    }

    b[sequence.ref_genome_length] = -b_0;

    return b;
}
std::vector<bam_api::ReadIndex> qmcp::SequentialCostScalingNetworkSolver::obtain_sequence(
    const bam_api::AOSPairedReads& sequence,
    const operations_research::SimpleMinCostFlow& min_cost_flow) {
    auto reduced_reads = std::vector<bam_api::ReadIndex>();

    for (std::size_t read_id = 0; read_id < sequence.reads.size(); ++read_id) {
        if (min_cost_flow.Flow(read_id) > 0) {
            reduced_reads.push_back(read_id);
        }
    }

    return reduced_reads;
}

