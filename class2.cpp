#include "class2.h"
#include "fastgr.h"
#include <iostream>

void Class2::make_chain(graph_t *g, node_t from, node_t to, edge_t orig)
{
    int r;
    node_t u, v;
    edge_t e;

    fastgraph fg(g, propmap);

    u = from;
    // we don't have label edges
    assert(boost::get(*to_virt_map, orig) == edge_t());
    for (r = boost::get(*rank_map, from) + 1; r <= boost::get(*rank_map, to); ++r) {
        if (r < boost::get(*rank_map, to)) {
	        // we create plain vnode only because we don't have label edges
            v = fg.virtual_node(); // don't care about increasing widths
            graph_t::children_iterator gi, gi_end;
            boost::tie(gi, gi_end) = g->children();
            for (int i = 0; gi != gi_end && i < 
                propmap.components[static_cast<int>(from)] + 1; ++gi, ++i)
                    boost::add_vertex(v, *gi);
            // add vertex_index map for v   
            boost::put(*rank_map, v, r);
        } else
            v = to;
        e = fg.virtual_edge(u, v, orig);
        // virtual weight is not important right now
        u = v;
    }
    assert(boost::get(*to_virt_map, orig) != edge_t());

}

void Class2::merge_chain(graph_t *g, edge_t e, edge_t f, int flag)
{
    edge_t rep;
    int lastrank = std::max(boost::get(*rank_map, boost::source(e, *g)), boost::get(*rank_map, boost::target(e, *g)));

    assert(boost::get(*to_virt_map, e) == edge_t());
    boost::put(*to_virt_map, e, f);
    rep = f;
    do {
        /* interclust multi-edges are not counted now */
        if (flag)
            boost::put(*count_map, rep, boost::get(*count_map, rep) + boost::get(*count_map, e));
        boost::put(*xpenalty_map, rep, boost::get(*xpenalty_map, rep) + boost::get(*xpenalty_map, e));
        boost::put(*weight_map, rep, boost::get(*weight_map, rep) + boost::get(*weight_map, e));
        if (boost::get(*rank_map, boost::target(rep, *g)) == lastrank)
            break;
	    // don't care about incr_width(g, aghead(rep));
        rep = boost::get(*out_map, boost::target(rep, *g))[0];

    } while (rep != edge_t());
}

void Class2::class2(graph_t * g)
{
    node_it n, ni_end;
    edge_t prev;

    for (boost::tie(n, ni_end) = vertices(*g); n != ni_end; ++n) {
	    // we don't use fast nodes, because they are the same as boost::vertices(*g)
        prev = edge_t();
        std::vector<edge_t> &out_n = boost::get(*orig_out_map, *n);
        for (std::vector<edge_t>::iterator eo = out_n.begin(); eo != out_n.end(); ++eo) {

            /* already processed */
            if (boost::get(*to_virt_map, *eo) != edge_t()) {
                prev = *eo;
                continue;
            }

            /* edges involving sub-clusters of g */

            fastgraph fg(g, propmap);
            /* merge multi-edges */
            if (prev != edge_t() && boost::source(*eo, *g) == boost::source(prev, *g) && 
                    boost::target(*eo, *g) == boost::target(prev, *g)) {
                if (boost::get(*rank_map, boost::source(*eo, *g)) == boost::get(*rank_map, boost::target(*eo, *g))) {
                    fg.merge_oneway(*eo, prev);
                    // currently, don't need to install other edges, make sure to use vector instead of boost::out_edges() when doing so
                    continue;
                }
                // We don't have label edges, so code flow definitely falls into this, we assume that ports_eq(e, prev) and our graph is not concentrated
                merge_chain(g, *eo, boost::get(*to_virt_map, prev), true);
                // currently, don't need to install other edges, make sure to use vector instead of boost::out_edges() when doing so
                continue;
                /* parallel edges with different labels fall through here */
            }

            /* self edges */
            if (boost::source(*eo, *g) == boost::target(*eo, *g)) {
                // currently, don't need to install other edges, make sure to use vector instead of boost::out_edges() when doing so
                prev = *eo;
                continue;
            }
            // don't care about the UF_find stuff

            /* flat edges */
            if (boost::get(*rank_map, boost::source(*eo, *g)) == boost::get(*rank_map, boost::target(*eo, *g))) {
                fg.flat_edge(g, *eo);
                prev = *eo;
                continue;
            }

            /* forward edges */
            if (boost::get(*rank_map, boost::target(*eo, *g)) > boost::get(*rank_map, boost::source(*eo, *g))) {
                make_chain(g, boost::source(*eo, *g), boost::target(*eo, *g), *eo);
                prev = *eo;
                continue;
            }

    	    /* backward edges */
            // suppose we only have forward edges
        }
    }

}

