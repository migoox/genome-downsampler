#include "../include/qmcp-solver/sequential_max_flow_solver.hpp"

#include <ortools/graph/max_flow.h>

#include <cstdint>
#include <optional>
#include <vector>

#include "bam-api/bam_paired_reads.hpp"

void create_network_flow_graph(operations_research::SimpleMaxFlow& max_flow,
                               const bam_api::AOSPairedReads& sequence,
                               unsigned int M);
std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence,
                                        unsigned int M);
std::vector<bam_api::ReadIndex> obtain_sequence(
    const bam_api::AOSPairedReads& sequence,
    const operations_research::SimpleMaxFlow& max_flow);
void add_pairs(
    std::vector<bam_api::ReadIndex>& reduced_reads,
    const std::vector<bool>& mapped_reads,
    const std::vector<std::optional<bam_api::ReadIndex>>& read_pair_map);

void qmcp::SequentialMaxFlowSolver::solve() {
    operations_research::SimpleMaxFlow max_flow;

    create_network_flow_graph(max_flow, input_sequence_, M_);

    int status = max_flow.Solve(input_sequence_.ref_genome_length + 1,
                                input_sequence_.ref_genome_length + 2);

    if (status == operations_research::SimpleMaxFlow::OPTIMAL) {
        output_sequence_ = obtain_sequence(input_sequence_, max_flow);
        std::cout << "OBTAINED FLOW == " << max_flow.OptimalFlow() << std::endl;
    } else {
        LOG(INFO) << "Solving the min cost flow problem failed. Solver status: "
                  << status;
    }
}

void create_network_flow_graph(operations_research::SimpleMaxFlow& max_flow,
                               const bam_api::AOSPairedReads& sequence,
                               unsigned int M) {
    std::vector<int64_t> start_nodes = {};
    std::vector<int64_t> end_nodes = {};
    std::vector<int64_t> capacities = {};

    // create normal edges
    for (int i = 0; i < sequence.reads.size(); ++i) {
        max_flow.AddArcWithCapacity(sequence.reads[i].start_ind - 1,
                                    sequence.reads[i].end_ind, 1);
    }

    // create backwards edges to push more flow
    for (int i = 0; i < sequence.ref_genome_length; ++i) {
        max_flow.AddArcWithCapacity(i + 1, i, INT64_MAX);
    }

    // create edge between our virtual node and 0
    // (0 -> sequence.ref_genome_length)
    max_flow.AddArcWithCapacity(0, sequence.ref_genome_length, INT64_MAX);

    // add supply and demand (negative supply = demand)
    std::vector<int> demand = create_demand_function(sequence, M);

    // create virtual sink and term
    uint64_t s = sequence.ref_genome_length + 1;
    uint64_t t = s + 1;

    for (int i = 0; i <= sequence.ref_genome_length; ++i) {
        if (demand[i] > 0)
            max_flow.AddArcWithCapacity(i, t, demand[i]);
        else if (demand[i] < 0)
            max_flow.AddArcWithCapacity(s, i, -demand[i]);
    }
}

std::vector<int> create_b_function(const bam_api::AOSPairedReads& sequence,
                                   unsigned int M) {
    std::vector<int> b(sequence.ref_genome_length + 1, 0);

    for (unsigned int i = 0; i < sequence.reads.size(); ++i) {
        for (unsigned int j = sequence.reads[i].start_ind;
             j < sequence.reads[i].end_ind; ++j) {
            b[j]++;
        }
    }

    // cap nucleotides with more reads than M to M
    for (unsigned int i = 0; i < sequence.ref_genome_length + 1; ++i) {
        if (b[i] > M) b[i] = M;
    }
    return b;
}

std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence,
                                        unsigned int M) {
    std::vector<int> b = create_b_function(sequence, M);

    int b_0 = b[0];
    for (int i = 0; i < sequence.ref_genome_length - 1; ++i) {
        b[i] = b[i] - b[i + 1];
    }

    b[sequence.ref_genome_length] = -b_0;

    return b;
}

std::vector<bam_api::ReadIndex> obtain_sequence(
    const bam_api::AOSPairedReads& sequence,
    const operations_research::SimpleMaxFlow& max_flow) {
    auto reduced_reads = std::vector<bam_api::ReadIndex>();
    std::vector<bool> mapped_reads =
        std::vector<bool>(sequence.read_pair_map.size());

    for (std::size_t read_id = 0; read_id < sequence.reads.size(); ++read_id) {
        if (max_flow.Flow(read_id) > 0) {
            reduced_reads.push_back(read_id);
            mapped_reads[read_id] = true;
        }
    }

    add_pairs(reduced_reads, mapped_reads, sequence.read_pair_map);
    return reduced_reads;
}

void add_pairs(
    std::vector<bam_api::ReadIndex>& reduced_reads,
    const std::vector<bool>& mapped_reads,
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

std::vector<bam_api::ReadIndex>
qmcp::SequentialMaxFlowSolver::output_sequence() {
    return output_sequence_;
}