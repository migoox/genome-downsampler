#include "../include/qmcp-solver/sequence_network_solver.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/named_function_params.hpp>
#include <iostream>
#include <boost/graph/cycle_canceling.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>

typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                           boost::property<boost::vertex_name_t, std::string>,
                           boost::property<boost::edge_capacity_t, long,
                           boost::property<boost::edge_residual_capacity_t, long,
                           boost::property<boost::edge_reverse_t, Traits::edge_descriptor>>> > Graph;

    // Define the properties for the edges
    typedef boost::property_map<Graph, boost::edge_capacity_t>::type CapacityMap;
    typedef boost::property_map<Graph, boost::edge_residual_capacity_t>::type ResidualCapacityMap;
    typedef boost::property_map<Graph, boost::edge_reverse_t>::type ReverseEdgeMap;

Graph create_circulation_Graph(const bam_api::BamSequence& sequence);

void qmcp::SequenceNetworkSolver::solve() {
    std::cout << "Not implemented!";
    // TODO(implement the function):
    // Boost graphs test:
    
    //typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS>
    //    Graph;

    Graph test_graph;
    boost::add_vertex(test_graph);
    boost::add_vertex(test_graph);
    boost::add_vertex(test_graph);

    boost::add_edge(0, 1, test_graph);
    boost::add_edge(1, 2, test_graph);



    //IMPLEMNATATION
    Graph network_graph = create_circulation_Graph(this->sequence);


    
}
Graph create_circulation_Graph(const bam_api::BamSequence& sequence) {
    Graph circulation(sequence.length);

    // TODO(borys): change it to be done in batch if possible
    boost::add_vertex(circulation);
    for(unsigned int i = 0; i< sequence.length; i++) {
        boost::add_edge(i, i+1, circulation);
    }

    // add edges
    for(auto read : sequence.reads) {
        boost::add_edge(read.start - 1, read.end, circulation);
    }


    boost::cycle_canceling(circulation,,,)




    boost::add_edge(0, 1, circulation);
    boost::add_edge(1, 2, circulation);
    return circulation;
}