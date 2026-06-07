#include "../include/qmcp-solver/quasi_mcp_cpu_max_flow_solver.hpp"

#include <cstddef>
#include <memory>
#include <vector>

#include "bam-api/region_api.hpp"
#include "logging/log.hpp"
#include "qmcp-solver/ortools_flow_types.hpp"
#include "qmcp-solver/solver.hpp"

std::unique_ptr<qmcp::Solution> qmcp::QuasiMcpCpuMaxFlowSolver::solve(uint32_t max_coverage,
                                                                      bam_api::RegionApi& reg_api) {
    input_sequence_ = reg_api.get_paired_reads_aos();

    operations_research::SimpleMaxFlow max_flow;

    create_network_flow_graph(max_flow, input_sequence_, max_coverage);

    const OrNode source = to_node(input_sequence_.ref_genome_length + 1);
    const OrNode sink = source + 1;
    const int status = max_flow.Solve(source, sink);

    if (status != operations_research::SimpleMaxFlow::OPTIMAL) {
        LOG_WITH_LEVEL(logging::INFO)
            << "Solving the min cost flow problem failed. Solver status: " << status;
    }

    return obtain_sequence(input_sequence_, max_flow);
}

void qmcp::QuasiMcpCpuMaxFlowSolver::create_network_flow_graph(
    operations_research::SimpleMaxFlow& max_flow, const bam_api::AOSPairedReads& sequence,
    unsigned int M) {
    const OrNode ref_len = to_node(sequence.ref_genome_length);
    const OrNode source = ref_len + 1;
    const OrNode sink = source + 1;

    for (const bam_api::Read& read : sequence.reads) {
        max_flow.AddArcWithCapacity(to_node(read.start_ind), to_node(read.end_ind + 1), to_flow(1));
    }

    for (OrNode i = 0; i < ref_len; ++i) {
        max_flow.AddArcWithCapacity(i + 1, i, std::numeric_limits<OrFlow>::max());
    }

    const std::vector<int> demand = create_demand_function(sequence, M);

    for (OrNode i = 0; i <= ref_len; ++i) {
        std::size_t ui = static_cast<std::size_t>(i);
        if (demand[ui] > 0) {
            max_flow.AddArcWithCapacity(i, sink, to_flow(demand[ui]));
        } else if (demand[ui] < 0) {
            max_flow.AddArcWithCapacity(source, i, to_flow(-demand[ui]));
        }
    }
}

std::vector<int> qmcp::QuasiMcpCpuMaxFlowSolver::create_b_function(
    const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int> b(sequence.ref_genome_length + 1, 0);

    for (std::size_t i = 0; i < sequence.reads.size(); ++i) {
        for (std::size_t j = sequence.reads[i].start_ind; j <= sequence.reads[i].end_ind; ++j) {
            ++b[j + 1];
        }
    }

    for (std::size_t i = 0; i < sequence.ref_genome_length + 1; ++i) {
        if (b[i] > static_cast<int>(M)) b[i] = static_cast<int>(M);
    }
    return b;
}

std::vector<int> qmcp::QuasiMcpCpuMaxFlowSolver::create_demand_function(
    const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int> b = create_b_function(sequence, M);

    const int b_1 = b[1];
    for (std::size_t i = 1; i < sequence.ref_genome_length; ++i) {
        b[i] = b[i] - b[i + 1];
    }

    b[0] = -b_1;

    return b;
}

std::unique_ptr<qmcp::Solution> qmcp::QuasiMcpCpuMaxFlowSolver::obtain_sequence(
    const bam_api::AOSPairedReads& sequence, const operations_research::SimpleMaxFlow& max_flow) {
    auto reduced_reads = std::make_unique<Solution>();

    for (bam_api::ReadIndex read_index = 0; read_index < sequence.reads.size(); ++read_index) {
        if (max_flow.Flow(to_arc(read_index)) > 0) {
            reduced_reads->push_back(read_index);
        }
    }

    return reduced_reads;
}
