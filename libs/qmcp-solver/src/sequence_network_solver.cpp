#include "../include/qmcp-solver/sequence_network_solver.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <iostream>

void qmcp::SequenceNetworkSolver::solve() {
    std::cout << "Not implemented!";
    // TODO(implement the function):
    // Boost graphs test:
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS>
        Graph;

    Graph test_graph;
    boost::add_vertex(test_graph);
    boost::add_vertex(test_graph);
    boost::add_vertex(test_graph);

    boost::add_edge(0, 1, test_graph);
    boost::add_edge(1, 2, test_graph);
}