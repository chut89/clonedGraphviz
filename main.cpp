#include "dotinit.h"
#include <iostream>
#include <cstdio>
//#include "types_bgl.h"

int main(int argc, char** argv) {

	graph_t g;
	typedef std::pair <int, int> E;
	E edge_array[] = 
//*/
{E(0, 1), E(2, 8), E(2, 9), E(3, 0), E(4, 0), E(5, 4), E(6, 3), E(7, 12), E(8, 19), E(9, 5), E(10, 7), E(10, 17), E(11, 13), E(12, 2), E(13, 6), E(13, 18), E(15, 14), E(15, 16), E(16, 1), E(17, 11), E(17, 11), E(17, 11), E(18, 5), E(20, 21), E(21, 22), E(22, 23)};//*/
   /*/
{E(0, 1), E(0, 2), E(0, 2), E(3, 4), E(5, 4)};//*/
    const int n_edges = sizeof(edge_array) / sizeof(E);
    std::vector<node_t> verts;
    for (std::size_t i = 0; i < 24; ++i) {
        verts.push_back(boost::add_vertex(g));
        //boost::put(boost::vertex_index, g, verts[i], i); 
    }
    for (std::size_t j = 0; j < n_edges; ++j) {
        bool inserted; edge_t tmp;
        boost::tie(tmp, inserted) = boost::add_edge(verts[edge_array[j].first], verts[edge_array[j].second], g); 
        //if (inserted) boost::put(boost::edge_index, g, tmp, j); 
    }

    dotinit dot;
    dot.dot_layout(&g);

    return 1;
} 



