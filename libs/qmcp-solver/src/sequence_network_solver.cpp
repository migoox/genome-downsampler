#include "../include/qmcp-solver/sequence_network_solver.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/properties.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <climits>
#include <iostream>
#include <boost/graph/cycle_canceling.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <vector>
#include "qmcp-solver/network_graph.hpp"


boost::NetworkGraph create_circulation_Graph(const bam_api::AOSPairedReads& sequence, unsigned int M);
std::vector<int> create_b_function(const bam_api::AOSPairedReads& sequence, unsigned int M);
std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence, unsigned int M);
bam_api::AOSPairedReads obtain_bamsequence(boost::Network::Graph& G,bam_api::AOSPairedReads & paired_reads);

void qmcp::SequenceNetworkSolver::solve() {
    std::cout << "Not implemented!";
    
    boost::NetworkGraph network_graph = create_circulation_Graph(this->sequence_, this->M_);

    boost::edmonds_karp_max_flow(network_graph.G, network_graph.s, network_graph.t);
    boost::cycle_canceling(network_graph.G);

    bam_api::AOSPairedReads output_sequence = obtain_bamsequence(network_graph.G,this->sequence_);
    
}

boost::NetworkGraph  create_circulation_Graph(const bam_api::AOSPairedReads& sequence, unsigned int M) {

    boost::Network::vertex_descriptor s;
    boost::Network::vertex_descriptor t;
    boost::Network::Graph g;

    boost::Network::size_type n(sequence.ref_genome_length + 1);

    for(boost::Network::size_type i = 0; i < n; ++i)   {
            add_vertex(g);
    }

    boost::Network::Capacity capacity = get(boost::edge_capacity, g);
    boost::Network::Reversed rev = get(boost::edge_reverse, g);
    boost::Network::ResidualCapacity residual_capacity = get(boost::edge_residual_capacity, g); 
    boost::Network::Weight weight = get(boost::edge_weight, g);
    boost::Network::Name name = get(boost::edge_name,g);

    boost::Network::EdgeAdder edge_adder(g, weight, capacity, rev, residual_capacity,name);

    // create backwards edges with infinity capacity
    for(unsigned int i = 0; i< sequence.ref_genome_length; i++) {
        edge_adder.addEdge(i + 1, i, 0, INT_MAX,-1);
    }

    //create normal edges
    for(unsigned int i = 0; i < sequence.reads.size(); i++) {
        int weight = (sequence.reads[i].start_ind - sequence.reads[i].end_ind) * sequence.reads[i].quality;
        edge_adder.addEdge(sequence.reads[i].start_ind - 1, sequence.reads[i].end_ind, weight, 1,i);
    }

    //create demand funciton
    std::vector<int> d = create_demand_function(sequence,M);


    //add s and t vertex
    // created edges from s to GRAPH and from t to GRAPH with correct capacites and weights
    s = add_vertex(g);
    t = add_vertex(g);

    for(unsigned int i = 1; i< sequence.ref_genome_length - 1; i++) {
        if(d[i] > 0) edge_adder.addEdge(i, t, 0, d[i],i);
        else if(d[i] < 0) edge_adder.addEdge(s, i, 0, -d[i],i); 
    }
    
    return boost::NetworkGraph{g,s,t};
}


std::vector<int> create_b_function(const bam_api::AOSPairedReads& sequence,unsigned int M) {
    
    std::vector<int> b(sequence.ref_genome_length, 0);

    for(unsigned int i =0; i<sequence.reads.size(); i++) {
        
        for(unsigned int j = sequence.reads[i].start_ind; j<sequence.reads[i].end_ind; j++) {
            b[j]++;
        }
    }

    //cap over M to M
    for(unsigned int i = 0; i<sequence.ref_genome_length; i++) {
        if(b[i] > M) b[i] = M;
    }

    return b;
}

std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence, unsigned int M) {
    std::vector<int> b = create_b_function(sequence,M);
    std::vector<int> d(sequence.ref_genome_length);

    for(int i =1; i<sequence.ref_genome_length - 1;i++) {
        d[i] = b[i] - b[i + 1];
    }
    d[0] = 0;
    return d;
}

bam_api::AOSPairedReads obtain_bamsequence(boost::Network::Graph& g,bam_api::AOSPairedReads & paired_reads){

    //std::vector<bam_api::Read> reads;
    bam_api::AOSPairedReads reduced_paired_reads;
    reduced_paired_reads.reads = std::vector<bam_api::Read>();;
    boost::Network::ResidualCapacity residual_capacity = get(boost::edge_residual_capacity, g);
    boost::Network::Weight weight = get(boost::edge_weight, g);
    boost::Network::Capacity capacity = get(boost::edge_capacity, g);
    boost::Network::Name name = get(boost::edge_name,g);

    for(auto edge_it = boost::edges(g).first; edge_it != boost::edges(g).second; ++edge_it) {
        auto edge = *edge_it;
        
        if(weight[edge] != 0 && residual_capacity[edge] < capacity[edge]) {
            bam_api::Read read = bam_api::Read {paired_reads.reads[name[edge]]};
            reduced_paired_reads.push_back(read);
        }
    }
    return reduced_paired_reads;
}   