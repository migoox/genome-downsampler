#include "../include/qmcp-solver/quasi_mcp_cpu_max_flow_solver.hpp"

#include <cstdint>
#include <memory>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "qmcp-solver/solver.hpp"

std::unique_ptr<qmcp::Solution> qmcp::QuasiMcpCpuMaxFlowSolver::solve(uint32_t max_coverage,
                                                                      bam_api::BamApi& bam_api) {
    input_sequence_ = bam_api.get_paired_reads_aos();

    operations_research::SimpleMaxFlow max_flow;

    create_network_flow_graph(max_flow, input_sequence_, max_coverage);

    int status = max_flow.Solve(input_sequence_.ref_genome_length + 1,
                                input_sequence_.ref_genome_length + 2);

    if (status != operations_research::SimpleMaxFlow::OPTIMAL) {
        LOG(INFO) << "Solving the min cost flow problem failed. Solver status: " << status;
    }

    return obtain_sequence(input_sequence_, max_flow);
}

void qmcp::QuasiMcpCpuMaxFlowSolver::create_network_flow_graph(
    operations_research::SimpleMaxFlow& max_flow, const bam_api::AOSPairedReads& sequence,
    unsigned int M) {
    // create normal edges
    for (const bam_api::Read& read : sequence.reads) {
        max_flow.AddArcWithCapacity(read.start_ind, read.end_ind + 1, 1);
    }

    // create backwards edges to push more flow
    for (int i = 0; i < sequence.ref_genome_length; ++i) {
        max_flow.AddArcWithCapacity(i + 1, i, INT64_MAX);
    }

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

std::vector<int> qmcp::QuasiMcpCpuMaxFlowSolver::create_b_function(
    const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int> b(sequence.ref_genome_length + 1, 0);

    for (unsigned int i = 0; i < sequence.reads.size(); ++i) {
        for (unsigned int j = sequence.reads[i].start_ind; j <= sequence.reads[i].end_ind; ++j) {
            ++b[j + 1];
        }
    }

    // cap nucleotides with more reads than M to M
    for (unsigned int i = 0; i < sequence.ref_genome_length + 1; ++i) {
        if (b[i] > M) b[i] = M;
    }
    return b;
}

std::vector<int> qmcp::QuasiMcpCpuMaxFlowSolver::create_demand_function(
    const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int> b = create_b_function(sequence, M);

    int b_1 = b[1];
    for (int i = 1; i < sequence.ref_genome_length; ++i) {
        b[i] = b[i] - b[i + 1];
    }

    b[0] = -b_1;

    return b;
}

std::unique_ptr<qmcp::Solution> qmcp::QuasiMcpCpuMaxFlowSolver::obtain_sequence(
    const bam_api::AOSPairedReads& sequence, const operations_research::SimpleMaxFlow& max_flow) {
    auto reduced_reads = std::make_unique<Solution>();

    for (bam_api::ReadIndex read_index = 0; read_index < sequence.reads.size(); ++read_index) {
        if (max_flow.Flow(read_index) > 0) {
            reduced_reads->push_back(read_index);
        }
    }

    return reduced_reads;
}
