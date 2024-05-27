#include "../include/qmcp-solver/sequential_max_flow_solver.hpp"

#include <cstdint>
#include <optional>
#include <vector>

#include "bam-api/bam_paired_reads.hpp"

qmcp::SequentialMaxFlowSolver::SequentialMaxFlowSolver() : is_data_loaded_(false){};

qmcp::SequentialMaxFlowSolver::SequentialMaxFlowSolver(const std::filesystem::path& filepath,
                                                       uint32_t min_seq_length,
                                                       uint32_t min_seq_mapq)
    : input_filepath_(filepath) {
    import_reads(input_filepath_, min_seq_length, min_seq_mapq);
}

void qmcp::SequentialMaxFlowSolver::solve(uint32_t max_coverage) {
    if (!is_data_loaded_) {
        std::cerr << "Couldn't run solver: data has not been loaded.\n";
        std::exit(EXIT_FAILURE);
    }
    operations_research::SimpleMaxFlow max_flow;

    create_network_flow_graph(max_flow, input_sequence_, max_coverage);

    int status = max_flow.Solve(input_sequence_.ref_genome_length + 1,
                                input_sequence_.ref_genome_length + 2);

    if (status == operations_research::SimpleMaxFlow::OPTIMAL) {
        output_sequence_ = obtain_sequence(input_sequence_, max_flow, find_pairs_);
    } else {
        LOG(INFO) << "Solving the min cost flow problem failed. Solver status: " << status;
    }
}

void qmcp::SequentialMaxFlowSolver::create_network_flow_graph(
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

std::vector<int> qmcp::SequentialMaxFlowSolver::create_b_function(
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

std::vector<int> qmcp::SequentialMaxFlowSolver::create_demand_function(
    const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int> b = create_b_function(sequence, M);

    int b_1 = b[1];
    for (int i = 1; i < sequence.ref_genome_length; ++i) {
        b[i] = b[i] - b[i + 1];
    }

    b[0] = -b_1;

    return b;
}

std::vector<bam_api::BAMReadId> qmcp::SequentialMaxFlowSolver::obtain_sequence(
    const bam_api::AOSPairedReads& sequence, const operations_research::SimpleMaxFlow& max_flow,
    bool find_pairs) {
    auto reduced_reads = std::vector<bam_api::BAMReadId>();

    std::vector<bool> mapped_reads;
    if (find_pairs) {
        mapped_reads.resize(sequence.read_pair_map.size(), false);
    }

    for (bam_api::ReadIndex read_index = 0; read_index < sequence.reads.size(); ++read_index) {
        if (max_flow.Flow(read_index) > 0) {
            bam_api::BAMReadId read_id = sequence.reads[read_index].bam_id;
            reduced_reads.push_back(read_id);
            if (find_pairs) {
                mapped_reads[read_id] = true;
            }
        }
    }

    if (find_pairs) {
        add_pairs(reduced_reads, mapped_reads, sequence.read_pair_map);
    }

    return reduced_reads;
}

void qmcp::SequentialMaxFlowSolver::find_pairs(bool flag) { find_pairs_ = flag; }

void qmcp::SequentialMaxFlowSolver::add_pairs(
    std::vector<bam_api::ReadIndex>& reduced_reads, const std::vector<bool>& mapped_reads,
    const std::vector<std::optional<bam_api::BAMReadId>>& read_pair_map) {
    for (const bam_api::BAMReadId& read_id : reduced_reads) {
        auto paired_read = read_pair_map[read_id];
        if (!paired_read.has_value()) continue;

        int pair_id = paired_read.value();
        if (!mapped_reads[pair_id]) {
            reduced_reads.push_back(pair_id);
        }
    }
}

std::vector<bam_api::BAMReadId> qmcp::SequentialMaxFlowSolver::output_sequence() {
    return output_sequence_;
}

void qmcp::SequentialMaxFlowSolver::export_reads(const std::filesystem::path& filepath) {
    bam_api::BamApi::write_bam(input_filepath_, filepath, output_sequence_);
}

void qmcp::SequentialMaxFlowSolver::import_reads(const std::filesystem::path& filepath,
                                                 uint32_t min_seq_length, uint32_t min_seq_mapq) {
    input_filepath_ = filepath;
    input_sequence_ = bam_api::BamApi::read_bam_aos(input_filepath_, min_seq_length, min_seq_mapq);
    is_data_loaded_ = true;
}

void qmcp::SequentialMaxFlowSolver::set_reads(const bam_api::AOSPairedReads& input_sequence) {
    input_sequence_ = input_sequence;
    is_data_loaded_ = true;
}

const std::vector<bam_api::BAMReadId>& qmcp::SequentialMaxFlowSolver::get_output() {
    return output_sequence_;
}