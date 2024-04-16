#include "../include/qmcp-solver/sequence_network_solver.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/named_function_params.hpp>
#include <climits>
#include <iostream>
#include <boost/graph/cycle_canceling.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <vector>
#include "qmcp-solver/network_graph.hpp"


boost::NetworkGraph::Graph create_circulation_Graph(const bam_api::BamSequence& sequence);
std::vector<int> create_b_function(const bam_api::BamSequence& sequence);

void qmcp::SequenceNetworkSolver::solve() {
    std::cout << "Not implemented!";
    
    boost::NetworkGraph::Graph network_graph = create_circulation_Graph(this->sequence);
    
}

boost::NetworkGraph::Graph  create_circulation_Graph(const bam_api::BamSequence& sequence) {

    boost::NetworkGraph::vertex_descriptor s;
    boost::NetworkGraph::vertex_descriptor t;
    boost::NetworkGraph::Graph g;

    boost::NetworkGraph::size_type n(sequence.length + 1);

    for(boost::NetworkGraph::size_type i = 0; i < n; ++i)   {
            add_vertex(g);
    }

    boost::NetworkGraph::Capacity capacity = get(boost::edge_capacity, g);
    boost::NetworkGraph::Reversed rev = get(boost::edge_reverse, g);
    boost::NetworkGraph::ResidualCapacity residual_capacity = get(boost::edge_residual_capacity, g); 
    boost::NetworkGraph::Weight weight = get(boost::edge_weight, g);

    boost::NetworkGraph::EdgeAdder edge_adder(g, weight, capacity, rev, residual_capacity);

    // create backwards edges with infinity capacity
    for(unsigned int i = 0; i< sequence.length; i++) {
        edge_adder.addEdge(i + 1, i, 0, INT_MAX);
    }

    //create normal edges
    for(unsigned int i = 0; i < sequence.reads.size(); i++) {
        int weight = (sequence.reads[i].start - sequence.reads[i].end) * sequence.reads[i].quality;
        edge_adder.addEdge(sequence.reads[i].start - 1, sequence.reads[i].end, weight, 1);
    }

    //create b funciton
    std::vector<int> b = create_b_function(sequence);


    return g;
}


std::vector<int> create_b_function(const bam_api::BamSequence& sequence, int M) {
    
    std::vector<int> b(sequence.length, 0);

    for(unsigned int i =0; i<sequence.reads.size(); i++) {
        
        for(unsigned int j = sequence.reads[i].start; j<sequence.reads[i].end; j++) {
            b[j]++;
        }
    }

    //cap over M to M
    for(unsigned int i = 0; i<sequence.length; i++) {
        if(b[i] > M) b[i] = M;
    }

    return b;
}