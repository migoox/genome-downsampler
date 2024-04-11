#include "../include/qmcp-solver/sequence_network_solver.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/named_function_params.hpp>
#include <iostream>
#include <boost/graph/cycle_canceling.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <vector>

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

    typedef std::pair<int, int> Edge;

    std::vector<Edge> edges(sequence.length + sequence.reads.size());
    for(unsigned int i = 0; i< sequence.length; i++) {
        //boost::add_edge(i, i+1, circulation);
        edges[i] = Edge{i,i+1};
    }

    for(unsigned int i = 0; i < sequence.reads.size(); i++) {
        edges[i + sequence.length] = Edge{sequence.reads[i].start - 1, sequence.reads[i].end};
    }

    Graph circulation(edges.data(),edges.data() + edges.size() / sizeof(Edge),sequence.length + 1);

    boost::cycle_canceling(circulation,,,)

    return circulation;
}