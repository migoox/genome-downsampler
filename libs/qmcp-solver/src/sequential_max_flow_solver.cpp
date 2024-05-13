#include "../include/qmcp-solver/sequential_max_flow_solver.hpp"

#include <ortools/graph/max_flow.h>

#include <cstdint>
#include <vector>

#include "bam-api/bam_paired_reads.hpp"

void create_network_flow_graph(
    operations_research::SimpleMaxFlow& min_cost_flow,
    const bam_api::AOSPairedReads& sequence, unsigned int M);
std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence,
                                        unsigned int M);
bam_api::AOSPairedReads obtain_sequence(
    const bam_api::AOSPairedReads& sequence,
    const operations_research::SimpleMaxFlow& min_cost_flow);

void qmcp::SequentialMaxFlowSolver::solve() {
    // TODO(implement the function):
    // ortools max flow example basing on the link:
    // https://developers.google.com/optimization/flow/maxflow#c++
    // Worth reading:
    // https://github.com/google/or-tools/blob/stable/ortools/graph/min_cost_flow.h
    operations_research::SimpleMaxFlow min_cost_flow;

    create_network_flow_graph(min_cost_flow, input_sequence_, M_);

    int status = min_cost_flow.Solve(input_sequence_.ref_genome_length + 1,
                                     input_sequence_.ref_genome_length + 2);

    if (status == 0) {
        output_sequence_ = obtain_sequence(input_sequence_, min_cost_flow);
    } else {
        LOG(INFO) << "Solving the min cost flow problem failed. Solver status: "
                  << status;
    }
}

void create_network_flow_graph(
    operations_research::SimpleMaxFlow& min_cost_flow,
    const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int64_t> start_nodes = {};
    std::vector<int64_t> end_nodes = {};
    std::vector<int64_t> capacities = {};

    // create normal edges
    for (int i = 0; i < sequence.reads.size(); ++i) {
        int weight = (sequence.reads[i].start_ind - sequence.reads[i].end_ind) *
                     sequence.reads[i].quality;
        min_cost_flow.AddArcWithCapacity(sequence.reads[i].start_ind - 1,
                                         sequence.reads[i].end_ind, 1);
    }

    // create backwards edge
    for (int i = 0; i < sequence.ref_genome_length; ++i) {
        min_cost_flow.AddArcWithCapacity(i + 1, i, M * 100);
    }

    min_cost_flow.AddArcWithCapacity(0, sequence.ref_genome_length, M * 100);

    // add supply and demand (negative supply = demand)
    std::vector<int> demand = create_demand_function(sequence, M);

    uint64_t s = sequence.ref_genome_length + 1;
    uint64_t t = s + 1;
    for (int i = 0; i < sequence.ref_genome_length; ++i) {
        // min_cost_flow.SetNodeSupply(i, demand[i]);
        if (demand[i] > 0)
            min_cost_flow.AddArcWithCapacity(i, t, demand[i]);
        else
            min_cost_flow.AddArcWithCapacity(s, i, -demand[i]);
    }
}

std::vector<int> create_b_function(const bam_api::AOSPairedReads& sequence,
                                   unsigned int M) {
    std::vector<int> b(sequence.ref_genome_length + 1, 0);

    for (unsigned int i = 0; i < sequence.reads.size(); ++i) {
        for (unsigned int j = sequence.reads[i].start_ind;
             j < sequence.reads[i].end_ind; j++) {
            if (j > sequence.ref_genome_length) {
                int p = 0;
                continue;
            }
            b[j]++;
        }
    }

    // cap nucleotides with more reads than M to M
    for (unsigned int i = 0; i < sequence.ref_genome_length; ++i) {
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

bam_api::AOSPairedReads obtain_sequence(
    const bam_api::AOSPairedReads& sequence,
    const operations_research::SimpleMaxFlow& min_cost_flow) {
    bam_api::AOSPairedReads reduced_paired_reads;
    reduced_paired_reads.reads = std::vector<bam_api::Read>();
    std::vector<bool> mapped_reads =
        std::vector<bool>(sequence.read_pair_map.size());

    for (std::size_t i = 0; i < sequence.reads.size(); ++i) {
        if (min_cost_flow.Flow(i) > 0) {
            bam_api::Read read = bam_api::Read{sequence.reads[i]};
            reduced_paired_reads.push_back(read);
            mapped_reads[read.id] = true;
        }
    }

    for (const bam_api::Read& read : reduced_paired_reads.reads) {
        auto paired_read = sequence.read_pair_map[read.id];
        if (paired_read.has_value()) {
            int pair_id = paired_read.value();
            if (!mapped_reads[paired_read.value()]) {
                mapped_reads[paired_read.value()] = true;
                reduced_paired_reads.push_back(
                    sequence.reads[paired_read.value()]);
            }
        }
    }

    return reduced_paired_reads;
}