#include "rank.h"
#include "ns.h"
#include "fastgr.h"
#include "decomp.h"
#include "class1.h"

#include <iostream>
#include <limits>

void Rank::dot_rank(graph_t *g) {

    Class1 cl(g, propmap);
    cl.class1(g);
	//acyclic(g);
	decomp::decompose(g, propmap.components);
	rank1(g);
    expand_ranksets(g);
    cleanup1(g);
}

void Rank::rank1(graph_t *g) {

    boost::ref_property_map<graph_t*, Agraphinfo_t>
        rootg_propt(boost::get_property(*g, graph_IDproperty));

	graph_t::children_iterator gi, gi_end;
    for (boost::tie(gi, gi_end) = g->children(); gi != gi_end; ++gi) {
        NetworkSimplex ns(&(*gi), propmap);
		//ns.dump_graph(&(*(gi))); 
		ns.rank(&(*gi), 1, 10);	

    }

}

void Rank::expand_ranksets(graph_t * g)
{
    int c;
    node_it n, ni_end;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
        graph_propt(boost::get_property(*g, graph_IDproperty));

    graph_propt[g].minrank = std::numeric_limits<short>::max(); // ugly here, when number of nodes > 2^15
    graph_propt[g].maxrank = -1;

    for (boost::tie(n, ni_end) = boost::vertices(*g); n != ni_end && *n != boost::graph_traits<graph_t>::null_vertex(); ++n) {
        // we don't have leader nodes, so don't care about these things
        const int &rank_n = boost::get(*rank_map, *n);
        if (graph_propt[g].maxrank < rank_n)
            graph_propt[g].maxrank = rank_n;
        if (graph_propt[g].minrank > rank_n)
            graph_propt[g].minrank = rank_n;

	}
    // we don't have clusters, so don't worry about set_minmax for them
}

void Rank::cleanup1(graph_t * g)
{
    node_it n, vi_end;
    out_edge_iterator eo, eo_nxt, eo_end;
    int c;

    mark_map = boost::get(&Agnodeinfo_t::mark, *g);

    for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) {
        boost::tie(eo, eo_end) = boost::out_edges(*n, *g);
        for (eo_nxt = eo; eo != eo_end; eo = eo_nxt) {
            ++eo_nxt;
            if (boost::get(*edge_type_map, *eo) == NORMAL)
                boost::put(*to_virt_map, *eo, edge_t());
            else if (boost::get(*edge_type_map, *eo) == VIRTUAL);
                boost::remove_edge(*eo, *g); // notice that edge iterators/descriptors invalid here
        }
    }

    for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) {
        boost::get(*in_map, *n).clear();
        boost::get(*out_map, *n).clear();
        boost::put(mark_map, *n, false);

    }


}
