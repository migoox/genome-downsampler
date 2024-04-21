
#ifndef NETWORK_GRAPH_DIRECTED_HPP
#define NETWORK_GRAPH_DIRECTED_HPP 

#include <boost/graph/properties.hpp>
#include <boost/pending/property.hpp>
#include <iostream>
#include <cstdlib>
#include <boost/graph/adjacency_list.hpp>

namespace boost {
    struct Network {
        typedef adjacency_list_traits < vecS, vecS, directedS > Traits;

        typedef adjacency_list < vecS, vecS, directedS, no_property,
                property < edge_capacity_t, long,
                    property < edge_residual_capacity_t, long,
                        property < edge_reverse_t, Traits::edge_descriptor, 
                            property <edge_weight_t, long,
                                property<edge_name_t, long>
                                > >
                            > 
                        > > Graph;
        typedef property_map < Graph, edge_capacity_t >::type Capacity;
        typedef property_map < Graph, edge_residual_capacity_t >::type ResidualCapacity;
        typedef property_map < Graph, edge_weight_t >::type Weight;
        typedef property_map < Graph, edge_reverse_t>::type Reversed;
        typedef property_map < Graph, edge_name_t>::type Name;
        typedef boost::graph_traits<Graph>::vertices_size_type size_type;
        typedef Traits::vertex_descriptor vertex_descriptor;
        
        class EdgeAdder {
        public:
            EdgeAdder(Graph & g, Weight & w, Capacity & c, Reversed & rev, ResidualCapacity & residualCapacity, Name & name) 
                : m_g(g), m_w(w), m_cap(c), m_resCap(residualCapacity), m_rev(rev), m_name(name) {}
            void addEdge(vertex_descriptor v, vertex_descriptor w, long weight, long capacity,long id) {
                Traits::edge_descriptor e,f;
                e = add(v, w, weight, capacity,id);
                f = add(w, v, -weight, 0, id);
                m_rev[e] = f; 
                m_rev[f] = e; 
            }
        private:
            Traits::edge_descriptor add(vertex_descriptor v, vertex_descriptor w, long weight, long capacity,long id) {
                bool b;
                Traits::edge_descriptor e;
                boost::tie(e, b) = add_edge(vertex(v, m_g), vertex(w, m_g), m_g);
                if (!b) {
                std::cerr << "Edge between " << v << " and " << w << " already exists." << std::endl;
                std::abort();
                }
                m_cap[e] = capacity;
                m_w[e] = weight;
                m_resCap[e] = weight;
                m_name[e] = id;
                return e;
            }
            Graph & m_g;
            Weight & m_w;
            Capacity & m_cap;
            ResidualCapacity & m_resCap;
            Reversed & m_rev;
            Name & m_name;
        };
    };

    struct NetworkGraph {
        Network::Graph G;
        Network::vertex_descriptor s,t;
    };
}
#endif
