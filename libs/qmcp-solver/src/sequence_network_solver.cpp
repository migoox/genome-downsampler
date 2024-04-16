#include "../include/qmcp-solver/sequence_network_solver.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/graph/named_function_params.hpp>
#include <iostream>
#include <boost/graph/cycle_canceling.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <vector>


    typedef boost::adjacency_list_traits <boost::vecS, boost::vecS, boost::directedS> Traits;

    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::directedS, boost::no_property,
            boost::property <boost::edge_capacity_t, long,
                boost::property <boost::edge_residual_capacity_t, long,
                    boost::property <boost::edge_reverse_t, Traits::edge_descriptor, 
                        boost::property <boost::edge_weight_t, long>
                             > 
                        > 
                     > > Graph;
    typedef boost::property_map < Graph, boost::edge_capacity_t >::type Capacity;
    typedef boost::property_map < Graph, boost::edge_residual_capacity_t >::type ResidualCapacity;
    typedef boost::property_map < Graph, boost::edge_weight_t >::type Weight;
    typedef boost::property_map < Graph, boost::edge_reverse_t>::type Reversed;
    typedef boost::graph_traits<Graph>::vertices_size_type size_type;
    typedef Traits::vertex_descriptor vertex_descriptor;

      class EdgeAdder {
    public:
        EdgeAdder(Graph & g, Weight & w, Capacity & c, Reversed & rev, ResidualCapacity & residualCapacity) 
            : m_g(g), m_w(w), m_cap(c), m_resCap(residualCapacity), m_rev(rev) {}
        void addEdge(vertex_descriptor v, vertex_descriptor w, long weight, long capacity) {
            Traits::edge_descriptor e,f;
            e = add(v, w, weight, capacity);
            f = add(w, v, -weight, 0);
            m_rev[e] = f; 
            m_rev[f] = e; 
        }
    private:
        Traits::edge_descriptor add(vertex_descriptor v, vertex_descriptor w, long weight, long capacity) {
            bool b;
            Traits::edge_descriptor e;
            boost::tie(e, b) = add_edge(vertex(v, m_g), vertex(w, m_g), m_g);
            if (!b) {
              std::cerr << "Edge between " << v << " and " << w << " already exists." << std::endl;
              std::abort();
            }
            m_cap[e] = capacity;
            m_w[e] = weight;
            return e;
        }
        Graph & m_g;
        Weight & m_w;
        Capacity & m_cap;
        ResidualCapacity & m_resCap;
        Reversed & m_rev;
    };


Graph create_circulation_Graph(const bam_api::BamSequence& sequence);

void qmcp::SequenceNetworkSolver::solve() {
    std::cout << "Not implemented!";
    

    
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