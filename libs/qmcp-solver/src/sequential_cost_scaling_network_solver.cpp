#include "../include/qmcp-solver/sequential_cost_scaling_network_solver.hpp"

#include <cstdint>
#include <iostream>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"

qmcp::SequentialCostScalingNetworkSolver::SequentialCostScalingNetworkSolver()
    : is_data_loaded_(false){};

qmcp::SequentialCostScalingNetworkSolver::SequentialCostScalingNetworkSolver(
    const std::filesystem::path& filepath, uint32_t min_seq_length, uint32_t min_seq_mapq)
    : input_filepath_(filepath) {
    import_reads(input_filepath_, min_seq_length, min_seq_mapq);
}

void qmcp::SequentialCostScalingNetworkSolver::solve(uint32_t max_coverage) {
    if (!is_data_loaded_) {
        std::cerr << "Couldn't run solver: data has not been loaded.\n";
        std::exit(EXIT_FAILURE);
    }
    operations_research::SimpleMinCostFlow min_cost_flow;

    create_network_flow_graph(min_cost_flow, input_sequence_, max_coverage);

    int status = min_cost_flow.Solve();

    if (status == operations_research::MinCostFlow::OPTIMAL) {
        output_sequence_ = obtain_sequence(input_sequence_, min_cost_flow);
    } else {
        LOG(INFO) << "Solving the min cost flow problem failed. Solver status: " << status;
    }
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
    std::vector<bool> mapped_reads = std::vector<bool>(sequence.read_pair_map.size());

    for (std::size_t read_id = 0; read_id < sequence.reads.size(); ++read_id) {
        if (min_cost_flow.Flow(read_id) > 0) {
            reduced_reads.push_back(read_id);
            mapped_reads[read_id] = true;
        }
    }

    add_pairs(reduced_reads, mapped_reads, sequence.read_pair_map);
    return reduced_reads;
}

void qmcp::SequentialCostScalingNetworkSolver::add_pairs(
    std::vector<bam_api::ReadIndex>& reduced_reads, const std::vector<bool>& mapped_reads,
    const std::vector<std::optional<bam_api::ReadIndex>>& read_pair_map) {
    for (const bam_api::ReadIndex& read_id : reduced_reads) {
        auto paired_read = read_pair_map[read_id];
        if (!paired_read.has_value()) continue;

        int pair_id = paired_read.value();
        if (!mapped_reads[pair_id]) {
            reduced_reads.push_back(pair_id);
        }
    }
}

std::vector<bam_api::ReadIndex> qmcp::SequentialCostScalingNetworkSolver::output_sequence() {
    return output_sequence_;
}

void qmcp::SequentialCostScalingNetworkSolver::export_reads(const std::filesystem::path& filepath) {
    bam_api::BamApi::write_bam(input_filepath_, filepath, output_sequence_);
}

void qmcp::SequentialCostScalingNetworkSolver::import_reads(const std::filesystem::path& filepath,
                                                            uint32_t min_seq_length,
                                                            uint32_t min_seq_mapq) {
    input_filepath_ = filepath;
    input_sequence_ = bam_api::BamApi::read_bam_aos(input_filepath_, min_seq_length, min_seq_mapq);
    is_data_loaded_ = true;
}
