#include "../include/qmcp-solver/sequential_cost_scaling_network_solver.hpp"

#include <cstdint>
#include <vector>


void create_network_flow_graph(operations_research::SimpleMinCostFlow& min_cost_flow, const bam_api::AOSPairedReads& sequence, unsigned int M);
std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence, unsigned int M);
bam_api::AOSPairedReads obtain_sequence(const bam_api::AOSPairedReads& sequence, const operations_research::SimpleMinCostFlow& min_cost_flow);

void qmcp::SequentialCostScalingNetworkSolver::solve() {
    // TODO(implement the function):
    // ortools max flow example basing on the link:
    // https://developers.google.com/optimization/flow/maxflow#c++
    // Worth reading:
    // https://github.com/google/or-tools/blob/stable/ortools/graph/min_cost_flow.h
    operations_research::SimpleMinCostFlow min_cost_flow;

    create_network_flow_graph(min_cost_flow, input_sequence_, M_);

    int status = min_cost_flow.Solve();

    if (status == operations_research::MinCostFlow::OPTIMAL) {
        output_sequence_ = obtain_sequence(input_sequence_, min_cost_flow);
    } else {
        LOG(INFO) << "Solving the min cost flow problem failed. Solver status: "
              << status;
    }
}


void create_network_flow_graph(operations_research::SimpleMinCostFlow& min_cost_flow, const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int64_t> start_nodes = {};
    std::vector<int64_t> end_nodes = {};
    std::vector<int64_t> capacities = {};

    //create normal edges
    for(int i = 0; i < sequence.reads.size(); i++) {
        int weight = (sequence.reads[i].start_ind - sequence.reads[i].end_ind) * sequence.reads[i].quality;
        min_cost_flow.AddArcWithCapacityAndUnitCost(sequence.reads[i].start_ind - 1, sequence.reads[i].end_ind, 1, weight);
    }

    //create backwards edge
    for(int i = 0; i< sequence.ref_genome_length; i++) {
        min_cost_flow.AddArcWithCapacityAndUnitCost(i + 1, i, M * 100, 0);
    }

    // add supply and demand (negative supply = demand)
    std::vector<int> demand = create_demand_function(sequence,M);
    for (int i = 0; i < sequence.ref_genome_length; ++i) {
        min_cost_flow.SetNodeSupply(i, demand[i]);
  }

}

std::vector<int> create_b_function(const bam_api::AOSPairedReads& sequence,unsigned int M) {
    std::vector<int> b(sequence.ref_genome_length, 0);

    for(unsigned int i =0; i<sequence.reads.size(); i++) {        
        for(unsigned int j = sequence.reads[i].start_ind; j<sequence.reads[i].end_ind; j++) {
            if(j > sequence.ref_genome_length) continue;
            b[j]++;
        }
    }

    //cap nucleotides with more reads than M to M
    for(unsigned int i = 0; i<sequence.ref_genome_length; i++) {
        if(b[i] > M) b[i] = M;
    }
    return b;
}

std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int> b = create_b_function(sequence,M);

    for(int i = 1; i<sequence.ref_genome_length - 1; i++) {
        b[i] = b[i] - b[i + 1];
    }
    b[0] = -b[0];

    return b;
}

bam_api::AOSPairedReads obtain_sequence(const bam_api::AOSPairedReads& sequence, const operations_research::SimpleMinCostFlow& min_cost_flow) {
    bam_api::AOSPairedReads reduced_paired_reads;
    reduced_paired_reads.reads = std::vector<bam_api::Read>();
    
    for (std::size_t i = 0; i < sequence.reads.size(); ++i) {
        int64_t cost = min_cost_flow.Flow(i) * min_cost_flow.UnitCost(i);
        if(min_cost_flow.Flow(i) > 0) {
            bam_api::Read read = bam_api::Read {sequence.reads[i]};
            reduced_paired_reads.push_back(read);
        }
    }
    
    return reduced_paired_reads;
}