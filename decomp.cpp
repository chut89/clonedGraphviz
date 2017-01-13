#include "decomp.h"

#include <boost/graph/copy.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include <iostream>

struct do_nothing
{
    template <typename VertexOrEdge1, typename VertexOrEdge2>
    void operator()(const VertexOrEdge1& , VertexOrEdge2& ) const {}
};

void decomp::decompose(graph_t *g, std::vector<int> &components) {
	undirected_graph_t copied_g;
	copy_graph(*g, copied_g, boost::vertex_copy(do_nothing()).edge_copy(do_nothing()));

    //std::vector<int> components(boost::num_vertices(copied_g));
    int num = boost::connected_components(copied_g, &components[0]);
	std::cout<<num<<std::endl;
	std::vector<graph_t*> g_comp = std::vector<graph_t*>(num);

	node_it n, n_end;
	int i = 0; 
 
	for (int i = 0; i < num; ++i)
		g_comp[i] = &((*g).create_subgraph());
	// the edges in subgraphs are induced automatically by the global graph
    for (boost::tie(n, n_end) = boost::vertices(*g); n != n_end; ++n, ++i) {
		boost::add_vertex(*n, *(g_comp[components[i]]));
        //boost::put(boost::vertex_index, *(g_comp[components[i]]), g_comp[components[i]]->global_to_local(*n), boost::num_vertices(*(g_comp[components[i]])) - 1); 
	}
}
