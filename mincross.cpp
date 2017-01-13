#include "mincross.h"
#include "class2.h"
#include "fastgr.h"
#include "decomp.h"
#include "timing.c"

#include <iostream>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <queue>

Mincross::Mincross(graph_t *g, PropertyMap propmap) : Root(g), propmap(propmap) {
    
    node_it n, vi_end;
    for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) {
        //boost::put(&Agnodeinfo_t::rank, *g, *n, boost::get(*(propmap.rank_map), *n));

        std::vector<edge_t> orig_in_n = boost::get(*(propmap.orig_in_map), *n);
        std::vector<edge_t> orig_in_n_loc;
        for (std::vector<edge_t>::iterator it = orig_in_n.begin(); it != orig_in_n.end(); ++it)
            orig_in_n_loc.push_back(*it);
        boost::put(&Agnodeinfo_t::orig_in, *g, *n, orig_in_n_loc);
        //*/
        std::vector<edge_t> orig_out_n = boost::get(*(propmap.orig_out_map), *n);
        std::vector<edge_t> orig_out_n_loc;
        for (std::vector<edge_t>::iterator it = orig_out_n.begin(); it != orig_out_n.end(); ++it)
            orig_out_n_loc.push_back(*it);
        boost::put(&Agnodeinfo_t::orig_out, *g, *n, orig_out_n_loc);//*/

        std::vector<edge_t> in_n = boost::get(*(propmap.in_map), *n);
        std::vector<edge_t> in_n_loc;
        for (std::vector<edge_t>::iterator it = in_n.begin(); it != in_n.end(); ++it)
            in_n_loc.push_back(*it);
        boost::put(&Agnodeinfo_t::in, *g, *n, in_n_loc);

        std::vector<edge_t> out_n = boost::get(*(propmap.out_map), *n);
        std::vector<edge_t> out_n_loc;
        for (std::vector<edge_t>::iterator it = out_n.begin(); it != out_n.end(); ++it)
            out_n_loc.push_back(*it);
        boost::put(&Agnodeinfo_t::out, *g, *n, out_n_loc);

        boost::put(&Agnodeinfo_t::node_type, *g, *n, NORMAL);
    }

    rank_map = boost::get(&Agnodeinfo_t::rank, *g);
    order_map = boost::get(&Agnodeinfo_t::order, *g);
    orig_in_map = boost::get(&Agnodeinfo_t::orig_in, *g);
    orig_out_map = boost::get(&Agnodeinfo_t::orig_out, *g);
    in_map = boost::get(&Agnodeinfo_t::in, *g);
    out_map = boost::get(&Agnodeinfo_t::out, *g);
    mark_map = boost::get(&Agnodeinfo_t::mark, *g);
    low_map = boost::get(&Agnodeinfo_t::low, *g);
    onstack_map = boost::get(&Agnodeinfo_t::onstack, *g);
    mval_map = boost::get(&Agnodeinfo_t::mval, *g);
    flat_in_map = boost::get(&Agnodeinfo_t::flat_in, *g);
    flat_out_map = boost::get(&Agnodeinfo_t::flat_out, *g);
    other_map = boost::get(&Agnodeinfo_t::other, *g);
    node_type_map = boost::get(&Agnodeinfo_t::node_type, *g);
    coord_map = boost::get(&Agnodeinfo_t::coord, *g);

    edge_it e, e_end;
    for (boost::tie(e, e_end) = boost::edges(*g); e != e_end; ++e) {
        boost::put(&Agedgeinfo_t::weight, *g, *e, boost::get(*(propmap.weight_map), *e));
        boost::put(&Agedgeinfo_t::xpenalty, *g, *e, boost::get(*(propmap.xpenalty_map), *e));
        boost::put(&Agedgeinfo_t::to_orig, *g, *e, boost::get(*(propmap.to_orig_map), *e));
        boost::put(&Agedgeinfo_t::to_virt, *g, *e, boost::get(*(propmap.to_virt_map), *e));
    }

    weight_map = boost::get(&Agedgeinfo_t::weight, *g);
    xpenalty_map = boost::get(&Agedgeinfo_t::xpenalty, *g);
    edge_type_map = boost::get(&Agedgeinfo_t::edge_type, *g);
    to_orig_map = boost::get(&Agedgeinfo_t::to_orig, *g);
    to_virt_map = boost::get(&Agedgeinfo_t::to_virt, *g);

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));
    graph_propt[g].rankdir = 1;
}

// g: global
void Mincross::mincross_options(graph_t * g)
{
    /* set default values */
    MinQuit = 8;
    MaxIter = 24;
    Convergence = .995;
}

/* allocate_ranks:
 * Allocate rank structure, determining number of nodes per rank.
 * Note that no nodes are put into the structure yet.
 */
// g: global
void Mincross::allocate_ranks(graph_t * g)
{
    int r, low, high, *cn;
    node_it n, vi_end;
    edge_t e;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    graph_propt[g].rank = new rank_t[graph_propt[g].maxrank + 2];

    for (r = 0; r < graph_propt[g].maxrank+2; ++r) 
        graph_propt[g].rank[r].flat = NULL;

    graph_t::children_iterator gi, gi_end;
    for (boost::tie(gi, gi_end) = g->children(); gi != gi_end; ++gi) {
        boost::ref_property_map<graph_t*, Agraphinfo_t>
                local_graph_propt(boost::get_property(*gi, graph_IDproperty));
        local_graph_propt[&(*gi)].rankdir = graph_propt[g].rankdir;
        // assume our graphs contain no flat edges
        local_graph_propt[&(*gi)].has_flat_edges = false;
        local_graph_propt[&(*gi)].rank = new rank_t[graph_propt[g].maxrank + 2];
        for (r = 0; r < graph_propt[g].maxrank+2; ++r) 
            local_graph_propt[&(*gi)].rank[r].flat = NULL;
        local_graph_propt[&(*gi)].minrank = graph_propt[g].minrank;
        local_graph_propt[&(*gi)].maxrank = graph_propt[g].maxrank;
    }

}

// g: global
void Mincross::init_mincross(graph_t * g)
{
    int size;

    if (Verbose)
        start_timer();

    ReMincross = false;
    Root = g;
    /* alloc +1 for the null terminator usage in do_ordering() */
    /* also, the +1 avoids attempts to alloc 0 sizes, something
       that efence complains about */
    node_it n, vi_end;
    for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) 
        size += boost::get(orig_out_map, *n).size();
    
    size++; // assume that g = dot_root(g)

    TE_list = new edge_t[size];
    TI_list = new int[size];
    mincross_options(g);
    // no need to do fillRanks(g);
    //*/
    Class2 cl(g, propmap);
    cl.class2(g);
    //decomp::decompose(g);//*/

    allocate_ranks(g);
    // we assume there is no graph and node ordering

}

// g: local
int Mincross::rcross(graph_t * g, int r)
{
    static int *Count, C;
    int top, bot, cross, max, i, k;
    node_t v;

    cross = 0;
    max = 0;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));
    //boost::ref_property_map<graph_t*, Agraphinfo_t>
    //        root_graph_propt(boost::get_property(*Root, graph_IDproperty));
    
    std::vector<node_t> &rtop = graph_propt[g].rank[r].v;

    int v_root_size = graph_propt[g].rank[r+1].v.size();
    if (C <= v_root_size) {
        C = v_root_size + 1;
        Count = new int[C];
    }

    for (i = 0; i < graph_propt[g].rank[r+1].v.size(); ++i)
        Count[i] = 0;

    for (std::vector<node_t>::iterator rtop_it = rtop.begin(); rtop_it != rtop.end(); ++rtop_it) {
        std::vector<edge_t> &rtop_out = boost::get(out_map, g->local_to_global(*rtop_it));
        if (max > 0) {
            for (std::vector<edge_t>::iterator it = rtop_out.begin(); it != rtop_out.end(); ++it) {
                for (k = boost::get(order_map, boost::target(*it, g->parent())) + 1; k <= max; ++k)
    		        cross += Count[k] * boost::get(xpenalty_map, *it);
            }
        }

        for (std::vector<edge_t>::iterator it = rtop_out.begin(); it != rtop_out.end(); ++it) {
            register int inv = boost::get(order_map, boost::target(*it, g->parent()));
            if (inv > max)
                max = inv;
            Count[inv] += boost::get(xpenalty_map, *it);
        }
    }

    // assume v has no port
    //std::vector<node_t> &rank_r_v = graph_propt[g].rank[r].v;
    //for (top = 0; top < rank_r_v.size(); ++top) {
    //    v = rank_r_v[top];
	    //if (ND_has_port(v))
	    //    cross += local_cross(ND_out(v), 1);
    //}
    // assume v has no port
    //std::vector<node_t> &rank_r_next_v = graph_propt[g].rank[r+1].v;
    //for (bot = 0; bot < rank_r_next_v.size(); bot++) {
        //v = rank_r_next_v[bot];
    	//if (ND_has_port(v))
        //    cross += local_cross(ND_in(v), -1);
    //}
    return cross;
}

// pass local graph to this method, global g and Root are just the same
int Mincross::ncross(graph_t * g)
{
    int r, count, nc;

    //g = Root;
    count = 0;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    for (r = graph_propt[g].minrank; r < graph_propt[g].maxrank; ++r) {
        if (graph_propt[g].rank[r].valid)
            count += graph_propt[g].rank[r].cache_nc;
	    else {
            nc = graph_propt[g].rank[r].cache_nc = rcross(g, r);
            count += nc;
            graph_propt[g].rank[r].valid = true;
        }
    }
    return count;
}
// g: local; v, w: local
void Mincross::exchange(graph_t *g, node_t v, node_t w)
{
    int vi, wi, r;

    r = boost::get(rank_map, g->local_to_global(v));
    vi = boost::get(order_map, g->local_to_global(v));
    wi = boost::get(order_map, g->local_to_global(w));

    boost::put(order_map, g->local_to_global(v), wi);

    //boost::ref_property_map<graph_t*, Agraphinfo_t>
    //        root_graph_propt(boost::get_property(*Root, graph_IDproperty));
    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    std::vector<node_t> &root_rank_r_v = graph_propt[g].rank[r].v;
    root_rank_r_v[wi] = v;

    boost::put(order_map, g->local_to_global(w), vi);
    
    root_rank_r_v[vi] = w;
}
// g: local, v: local, w: local
bool Mincross::left2right(graph_t * g, node_t v, node_t w)
{
    adjmatrix_t *M;
    bool rv;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(g->parent(), graph_IDproperty));
    /* CLUSTER indicates orig nodes of clusters, and vnodes of skeletons */
    if (!ReMincross) {
        // we assume our graph contains no cluster
    } else {
        // we assume our graph contains no cluster
    }
    M = graph_propt[g].rank[boost::get(rank_map, g->local_to_global(v))].flat;
    if (M == NULL)
        rv = false;
    else {
        if (graph_propt[g].rankdir & 1) {
	        node_t t = v;
	        v = w;
	        w = t;
        }
        int flat_idx_v = boost::get(low_map, g->local_to_global(v));
        int flat_idx_w = boost::get(low_map, g->local_to_global(w));
        rv = M->data[flat_idx_v * M->ncols + flat_idx_w];
    }
    return rv;
}
// v, w: global
int Mincross::in_cross(node_t v, node_t w)
{
    register int inv, cross = 0, t;

    std::vector<edge_t> &in_w = boost::get(in_map, w);
    for (std::vector<edge_t>::iterator it2 = in_w.begin(); it2 != in_w.end(); ++it2) {
        register int cnt = boost::get(xpenalty_map, *it2);

        inv = boost::get(order_map, boost::source(*it2, *Root));

        std::vector<edge_t> &in_v = boost::get(in_map, v);
        for (std::vector<edge_t>::iterator it1 = in_v.begin(); it1 != in_v.end(); ++it1) {
            t = boost::get(order_map, boost::source(*it1, *Root)) - inv;
            if (t > 0)
                cross += boost::get(xpenalty_map, *it1) * cnt;
        }
    }
    return cross;
}
// v, w: global
int Mincross::out_cross(node_t v, node_t w)
{
    register int inv, cross = 0, t;

    std::vector<edge_t> &out_w = boost::get(out_map, w);
    for (std::vector<edge_t>::iterator it2 = out_w.begin(); it2 != out_w.end(); ++it2) {
        register int cnt = boost::get(xpenalty_map, *it2);
        inv = boost::get(order_map, boost::target(*it2, *Root));

        std::vector<edge_t> &out_v = boost::get(out_map, v);
        for (std::vector<edge_t>::iterator it1 = out_v.begin(); it1 != out_v.end(); ++it1) {
            t = boost::get(order_map, boost::target(*it1, *Root)) - inv;
            if (t > 0)
                cross += boost::get(xpenalty_map, *it1) * cnt;
        }
    }
    return cross;

}

// g: local
int Mincross::transpose_step(graph_t * g, int r, bool reverse)
{
    int i, c0, c1, rv;
    node_t v, w;

    rv = 0;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    graph_propt[g].rank[r].candidate = false;

    std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
    std::vector<node_t>::iterator itStart = g_rank_r_v.begin(), itNext;

    itNext = itStart;
    ++itNext;
    for (; itNext != g_rank_r_v.end(); ++itStart, ++itNext) {
        v = *itStart;
        w = *itNext;
        assert(boost::get(order_map, g->local_to_global(v)) < boost::get(order_map, g->local_to_global(w)));
        if (left2right(g, v, w))
            continue;
        c0 = c1 = 0;
        if (r > 0) {
        	    c0 += in_cross(g->local_to_global(v), g->local_to_global(w));
        	    c1 += in_cross(g->local_to_global(w), g->local_to_global(v));
        }
        if (graph_propt[g].rank[r+1].v.size()) {
            c0 += out_cross(g->local_to_global(v), g->local_to_global(w));
            c1 += out_cross(g->local_to_global(w), g->local_to_global(v));
        }
        if ((c1 < c0) || ((c0 > 0) && reverse && (c1 == c0))) {
            exchange(g, v, w);
            rv += (c0 - c1);

            //boost::ref_property_map<graph_t*, Agraphinfo_t>
            //        root_graph_propt(boost::get_property(*g, graph_IDproperty));

            graph_propt[g].rank[r].valid = false;
            graph_propt[g].rank[r].candidate = true;

            if (r > graph_propt[g].minrank) {
                graph_propt[g].rank[r-1].valid = false;
                graph_propt[g].rank[r-1].candidate = true;
            }
            if (r < graph_propt[g].maxrank) {
                graph_propt[g].rank[r+1].valid = false;
                graph_propt[g].rank[r+1].candidate = true;
            }
        }
    }

    return rv;
}

// g: local
void Mincross::transpose(graph_t * g, bool reverse)
{
    int r, delta;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; ++r)
        graph_propt[g].rank[r].candidate = true;
    do {
        delta = 0;
#ifdef NOTDEF
	/* don't run both the upward and downward passes- they cancel. 
	   i tried making it depend on whether an odd or even pass, 
	   but that didn't help. 
	for (r = GD_maxrank(g); r >= GD_minrank(g); r--)
	    if (GD_rank(g)[r].candidate)
		delta += transpose_step(g, r, reverse); */
#endif
        for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; ++r) {
            if (graph_propt[g].rank[r].candidate)
                delta += transpose_step(g, r, reverse);
        }
    /*} while (delta > ncross(g)*(1.0 - Convergence)); */
    } while (delta >= 1);
}


// g: local
void Mincross::install_in_rank(graph_t * g, node_t n)
{
    int i, r;

    r = boost::get(rank_map, g->local_to_global(n));

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    i = graph_propt[g].rank[r].v.size();

    graph_propt[g].rank[r].v.push_back(n);
    boost::put(order_map, g->local_to_global(n), i);
#ifdef DEBUG1
    {
	node_it v, vi_end;

    for (boost::tie(v, vi_end) = boost::vertices(*g); v != vi_end; ++v)
        if (*v == n)
            break;
    assert(*v != boost::graph_traits<graph_t>::null_vertex());
    }
#endif

    //boost::ref_property_map<graph_t*, Agraphinfo_t>
    //        root_graph_propt(boost::get_property(*Root, graph_IDproperty));

}

// g: local
void Mincross::enqueue_neighbors(graph_t *g, std::queue<node_t> &q, node_t n0, int pass)
{
    int i;
    edge_t e;

    if (pass == 0) {
        std::vector<edge_t> &out_n0 = boost::get(out_map, g->local_to_global(n0));
        for (std::vector<edge_t>::iterator it = out_n0.begin(); it != out_n0.end(); ++it) {
            node_t head_e = boost::target(*it, g->parent());
            if (!boost::get(mark_map, head_e)) {
                boost::put(mark_map, head_e, true);
                q.push(g->global_to_local(head_e));
            }
        }
    } else {
        std::vector<edge_t> &in_n0 = boost::get(in_map, g->local_to_global(n0));
        for (std::vector<edge_t>::iterator it = in_n0.begin(); it != in_n0.end(); ++it) {
            node_t tail_e = boost::source(*it, g->parent());
            if (!boost::get(mark_map, tail_e)) {
                boost::put(mark_map, tail_e, true);
                q.push(g->global_to_local(tail_e));
            }
        }
    }
}

// g: local graph
void Mincross::build_ranks(graph_t * g, int pass) {
    int i, j;
    node_t n0;
    std::vector<edge_t> otheredges;
    std::queue<node_t> q;

    node_it n, vi_end;
    for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) 
    	boost::put(mark_map, g->local_to_global(*n), false);

#ifdef DEBUG1
    {
        edge_t e;
        for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) {    
            std::vector<edge_t> &out_n = boost::get(out_map, g->local_to_global(*n));
            for (std::vector<edge_t>::iterator it = out_n.begin(); it != out_n.end(); ++it)
                assert(boost::get(mark_map, boost::target(*it, g->parent())) == false);

            std::vector<edge_t> &in_n = boost::get(in_map, g->local_to_global(*n));
            for (std::vector<edge_t>::iterator it = in_n.begin(); it != in_n.end(); ++it)
                assert(boost::get(mark_map, boost::source(*it, g->parent())) == false);

        }
    }
#endif

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    for (i = graph_propt[g].minrank; i <= graph_propt[g].maxrank; ++i)
        graph_propt[g].rank[i].v.clear();
    std::cout<<"num_vertices in build_ranks:"<<boost::num_vertices(*g)<<std::endl;
    for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) {  
    	otheredges = ((pass == 0) ? boost::get(in_map, g->local_to_global(*n)) : 
                boost::get(out_map, g->local_to_global(*n)));
        if (otheredges.size())
            continue;
        if (boost::get(mark_map, g->local_to_global(*n)) == false) {
            boost::put(mark_map, g->local_to_global(*n), true);
            q.push(*n);
            while (!q.empty()) {
                node_t n0 = q.front();
                q.pop();
                // ranktype of n0 is not CLUSTER in our case
	            install_in_rank(g, n0);
	            enqueue_neighbors(g, q, n0, pass);
            }
        }
    }

    if (!q.empty())
        fprintf(stderr, "surprise\n");
    
    //boost::ref_property_map<graph_t*, Agraphinfo_t>
    //        root_graph_propt(boost::get_property(*Root, graph_IDproperty));
    for (i = graph_propt[g].minrank; i <= graph_propt[g].maxrank; ++i) {
        graph_propt[g].rank[i].valid = false;

        if ((graph_propt[g].rankdir & 1) && graph_propt[g].rank[i].v.size()) {
            int n, ndiv2;
            std::vector<node_t> &vlist = graph_propt[g].rank[i].v;
            n = vlist.size() - 1;
            ndiv2 = n / 2;
            for (j = 0; j <= ndiv2; j++)
                exchange(g, vlist[j], vlist[n - j]);
        }    
    }
    
    // in our case g == dot_root(g)
    if (ncross(g) > 0)
        transpose(g, false);
}
// g: local, v: local, w: local
void Mincross::balanceNodes(graph_t * g, int r, node_t v, node_t w)
{
    node_t s;			/* separator node */
    int sepIndex = 0;
    int nullType;		/* type of null nodes */
    int cntDummy = 0, cntOri = 0;
    int k = 0, m = 0, k1 = 0, m1 = 0, i = 0;

    /* we only consider v and w of different types */
    if (boost::get(node_type_map, g->local_to_global(v)) == boost::get(node_type_map, g->local_to_global(w)))
        return;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    /* count the number of dummy and original nodes */
    std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
    for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it) {
        if (boost::get(node_type_map, g->local_to_global(*it)) == NORMAL)
            cntOri++;
        else
            cntDummy++;
    }

    if (cntOri < cntDummy) {
        if (boost::get(node_type_map, g->local_to_global(v)) == NORMAL)
            s = v;
        else
            s = w;
    } else {
        if (boost::get(node_type_map, v) == NORMAL)
            s = w;
        else
            s = v;
    }

    /* get the separator node index */
    i = 0;
    for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it, ++i) {
        if (*it == s)
            sepIndex = i;
    }
    nullType = (boost::get(node_type_map, g->local_to_global(s)) == NORMAL) ? VIRTUAL : NORMAL;

    /* count the number of null nodes to the left and 
     * right of the separator node 
     */
    for (i = sepIndex - 1; i >= 0; i--) {
        if (boost::get(node_type_map, g->local_to_global(g_rank_r_v[i])) == nullType)
            k++;
        else
            break;
    }
    
    for (i = sepIndex + 1; i < g_rank_r_v.size(); i++) {
        if (boost::get(node_type_map, g->local_to_global(g_rank_r_v[i])) == nullType)
            m++;
        else
            break;
    }

    /* now exchange v,w and calculate the same counts */

    exchange(g, v, w);

    /* get the separator node index */
    i = 0;
    for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it, ++i) {
        if (*it == s)
            sepIndex = i;
    }

    /* count the number of null nodes to the left and 
     * right of the separator node 
     */
    for (i = sepIndex - 1; i >= 0; i--) {
        if (boost::get(node_type_map, g->local_to_global(g_rank_r_v[i])) == nullType)
            k1++;
        else
            break;
    }

    for (i = sepIndex + 1; i < g_rank_r_v.size(); i++) {
        if (boost::get(node_type_map, g->local_to_global(g_rank_r_v[i])) == nullType)
            m1++;
        else
            break;
    }

    if (abs(k1 - m1) > abs(k - m)) {
        exchange(g, v, w);		//revert to the original ordering
    }
}

// g: local
int Mincross::balance(graph_t * g)
{
    int i, c0, c1, rv;
    node_t v, w;
    int r;

    rv = 0;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(g->parent(), graph_IDproperty));

    for (r = graph_propt[g].maxrank; r >= graph_propt[g].minrank; --r) {

        graph_propt[g].rank[r].candidate = false;
        std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
        std::vector<node_t>::iterator itStart = g_rank_r_v.begin();
        std::vector<node_t>::iterator itNext = itStart;
        ++itNext;
        for (; itNext != g_rank_r_v.end(); ++itStart, ++itNext) {
            assert(boost::get(order_map, g->local_to_global(*itStart)) < boost::get(order_map, g->local_to_global(*itNext)));
            if (left2right(g, *itStart, *itNext))
                continue;
            c0 = c1 = 0;
            if (r > 0) {
                c0 += in_cross(g->local_to_global(*itStart), g->local_to_global(*itNext));
                c1 += in_cross(g->local_to_global(*itNext), g->local_to_global(*itStart));
            }
            std::vector<node_t> &g_rank_r_next_v = graph_propt[g].rank[r+1].v;
            if (g_rank_r_next_v.size() > 0) {
                c0 += out_cross(g->local_to_global(*itStart), g->local_to_global(*itNext));
                c1 += out_cross(g->local_to_global(*itNext), g->local_to_global(*itStart));
            }
        #if 0
            if ((c1 < c0) || ((c0 > 0) && reverse && (c1 == c0))) {
	        exchange(v, w);
	        rv += (c0 - c1);
	        GD_rank(Root)[r].valid = FALSE;
	        GD_rank(g)[r].candidate = TRUE;

	        if (r > GD_minrank(g)) {
	            GD_rank(Root)[r - 1].valid = FALSE;
	            GD_rank(g)[r - 1].candidate = TRUE;
	        }
	        if (r < GD_maxrank(g)) {
	            GD_rank(Root)[r + 1].valid = FALSE;
	            GD_rank(g)[r + 1].candidate = TRUE;
	        }
            }
        #endif

            if (c1 <= c0) {
                balanceNodes(g, r, *itStart, *itNext);
            }
        }
    }
    return rv;
}

// g: local
int Mincross::mincross(graph_t * g, int startpass, int endpass, int doBalance)
{
    int maxthispass, iter, trying, pass;
    int cur_cross, best_cross;

    //graph_t pg = g->parent();

    if (startpass > 1) {
        cur_cross = best_cross = ncross(g);
        save_best(g);
    } else
        cur_cross = best_cross = std::numeric_limits<short>::max();
    for (pass = startpass; pass <= endpass; pass++) {
        if (pass <= 1) {
            maxthispass = std::min(4, MaxIter);
            // always true in our case
            //if (g == dot_root(g))
            build_ranks(g, pass);

            if (pass == 0)
                flat_breakcycles(g);
            flat_reorder(g);

            if ((cur_cross = ncross(g)) <= best_cross) {
                save_best(g);
                best_cross = cur_cross;
            }
            trying = 0;
        } else {
            maxthispass = MaxIter;
            if (cur_cross > best_cross)
                restore_best(g);
            cur_cross = best_cross;
        }
        trying = 0;
        for (iter = 0; iter < maxthispass; iter++) {
            if (Verbose)
                fprintf(stderr,
                    "mincross: pass %d iter %d trying %d cur_cross %d best_cross %d\n",
                        pass, iter, trying, cur_cross, best_cross);
            if (trying++ >= MinQuit)
                break;
            if (cur_cross == 0)
                break;
            mincross_step(g, iter);
            if ((cur_cross = ncross(g)) <= best_cross) {
                save_best(g);
                if (cur_cross < Convergence * best_cross)
                    trying = 0;
                best_cross = cur_cross;
            }
        }
        if (cur_cross == 0)
            break;
    }
    if (cur_cross > best_cross) 
        restore_best(g);
    if (best_cross > 0) {
        transpose(g, false);
        best_cross = ncross(g);
    }
    if (doBalance) {
        for (iter = 0; iter < maxthispass; iter++)
            balance(g);
    }
//*/
    return best_cross;
}

struct less_than_pos {
    graph_t *g;
    less_than_pos(graph_t *g) : g(g) {}
    inline bool operator() (const node_t &n0, const node_t &n1)
    {
        order_map_t order_map = boost::get(&Agnodeinfo_t::order, g->parent());
        return boost::get(order_map, n0) - boost::get(order_map, n1);
    }

};

// g: local
void Mincross::restore_best(graph_t * g)
{
    node_t n;
    int i, r;

    /* for (n = GD_nlist(g); n; n = ND_next(n)) */
	/* ND_order(n) = saveorder(n); */
    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; ++r) {
        std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
        for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it)
            boost::put(order_map, g->local_to_global(*it), boost::get(coord_map, g->local_to_global(*it)).x);
    }

    //boost::ref_property_map<graph_t*, Agraphinfo_t>
    //        root_graph_propt(boost::get_property(*Root, graph_IDproperty));

    for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; ++r) {
        graph_propt[g].rank[r].valid = false;
        std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
        std::sort(g_rank_r_v.begin(), g_rank_r_v.end(), less_than_pos(g));

    }
}

// g: local
void Mincross::save_best(graph_t * g)
{
    node_t n;
    /* for (n = GD_nlist(g); n; n = ND_next(n)) */
	/* saveorder(n) = ND_order(n); */
    int i, r;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; ++r) {
        std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
        for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it)
            boost::get(coord_map, g->local_to_global(*it)).x = boost::get(order_map, g->local_to_global(*it));
    }
}
// g: local, n: local
bool Mincross::flat_mval(graph_t *g, node_t n)
{
    int i;
    node_t nn;

    if (boost::get(flat_in_map, g->local_to_global(n)).size()) {
        std::vector<edge_t> &fl = boost::get(flat_in_map, g->local_to_global(n));
        std::vector<edge_t>::iterator it = fl.begin();
        nn = boost::source(*it, g->parent());
        for (; it != fl.end(); ++it)
            if (boost::get(order_map, boost::source(*it, g->parent())) > boost::get(order_map, nn))
                nn = boost::source(*it, g->parent());
        if (boost::get(mval_map, nn) >= 0) {
            boost::put(mval_map, g->local_to_global(n), boost::get(mval_map, nn) + 1);
            return false;
        }
    } else if (boost::get(flat_out_map, g->local_to_global(n)).size()) {
        std::vector<edge_t> &fl = boost::get(flat_out_map, g->local_to_global(n));
        std::vector<edge_t>::iterator it = fl.begin();
        nn = boost::target(*it, g->parent());
        for (; it != fl.end(); ++it)
            if (boost::get(order_map, boost::target(*it, g->parent())) < boost::get(order_map, nn))
                nn = boost::target(*it, g->parent());
        if (boost::get(mval_map, nn) >= 0) {
            boost::put(mval_map, g->local_to_global(n), boost::get(mval_map, nn) - 1);
            return false;
        }
    }
    return true;
}

int ordercmpf(int *i0, int *i1)
{
    return (*i0) - (*i1);
}

// g: local
bool Mincross::medians(graph_t * g, int r0, int r1)
{

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    int i, j, j0, lm, rm, lspan, rspan, *list;
    bool hasfixed = false;

    list = TI_list;
    
    std::vector<node_t> &g_rank_r0_v = graph_propt[g].rank[r0].v;
    for (std::vector<node_t>::iterator it = g_rank_r0_v.begin(); it != g_rank_r0_v.end(); ++it) {
        //n = *it;
        j = 0;
        if (r1 > r0) {
            std::vector<edge_t> &out_n = boost::get(out_map, g->local_to_global(*it));
            for (std::vector<edge_t>::iterator n_it = out_n.begin(); n_it != out_n.end(); ++n_it) 
                // assume there is no port order
                if (boost::get(xpenalty_map, *n_it) > 0)
                    list[j++] = MC_SCALE * boost::get(order_map, boost::target(*n_it, g->parent()));
            
        } else {
            std::vector<edge_t> &in_n = boost::get(in_map, g->local_to_global(*it));
            for (std::vector<edge_t>::iterator n_it = in_n.begin(); n_it != in_n.end(); ++n_it) 
                // assume there is no port order
                if (boost::get(xpenalty_map, *n_it) > 0)
                    list[j++] = MC_SCALE * boost::get(order_map, boost::source(*n_it, g->parent()));
    
        }
        switch (j) {
            case 0:
                boost::put(mval_map, g->local_to_global(*it), -1);
                break;
            case 1:
                boost::put(mval_map, g->local_to_global(*it), list[0]);
                break;
            case 2:
                boost::put(mval_map, g->local_to_global(*it), (list[0] + list[1]) / 2);
                break;
            default:
                std::qsort(list, j, sizeof(int), (qsort_cmpf) ordercmpf);
                if (j % 2)
                    boost::put(mval_map, g->local_to_global(*it), list[j/2]);
                else {
                    /* weighted median */
                    rm = j / 2;
                    lm = rm - 1;
                    rspan = list[j - 1] - list[rm];
                    lspan = list[lm] - list[0];
                    if (lspan == rspan)
                        boost::put(mval_map, g->local_to_global(*it), (list[lm] + list[rm]) / 2);
                    else {
                        int w = list[lm] * rspan + list[rm] * lspan;
                        boost::put(mval_map, g->local_to_global(*it), w / (lspan + rspan));
                    }
                }
        }
    }

    for (std::vector<node_t>::iterator it = g_rank_r0_v.begin(); it != g_rank_r0_v.end(); ++it) {
        if (!boost::get(out_map, g->local_to_global(*it)).size() && !boost::get(in_map, g->local_to_global(*it)).size())
            hasfixed = hasfixed || flat_mval(g, *it);
    }
    return hasfixed;
}
// v: global
bool Mincross::is_a_normal_node_of(graph_t *g, node_t v) {
    return boost::get(node_type_map, v) == NORMAL;
}
// v: global
bool Mincross::is_a_vnode_of_an_edge_of(graph_t *g, node_t v) {
    if (boost::get(node_type_map, v) == VIRTUAL &&
        boost::get(in_map, v).size() == 1 && boost::get(out_map, v).size() == 1) {
    edge_t e = boost::get(out_map, v)[0];
    while (boost::get(edge_type_map, e) != NORMAL)
        e = boost::get(to_orig_map, e);
        return true;
    }
    return false;
}
// v: global
bool Mincross::inside_cluster(graph_t * g, node_t v)
{
    return (is_a_normal_node_of(g, v) | is_a_vnode_of_an_edge_of(g, v));
}

// g: local, v: local, e: global
bool Mincross::constraining_flat_edge(graph_t *g, node_t v, edge_t e)
{
    if (boost::get(weight_map, e) == 0) return false;    
	if (!inside_cluster(g, boost::target(e, g->parent()))) return false;
	if (!inside_cluster(g, boost::source(e, g->parent()))) return false;
	return true;
}

/* construct nodes reachable from 'here' in post-order.
* This is the same as doing a topological sort in reverse order.
*/
// g: local, v: local
int Mincross::postorder(graph_t * g, node_t v, node_t *list, int r)
{
    edge_t e;
    int i, cnt = 0;

    boost::put(mark_map, g->local_to_global(v), true);
    std::vector<edge_t> &flat_out_v = boost::get(flat_out_map, g->local_to_global(v));
    if (flat_out_v.size()) {
        for (std::vector<edge_t>::iterator it = flat_out_v.begin(); it != flat_out_v.end(); ++it) {
            if (!constraining_flat_edge(g, v, *it)) continue;
            if (!boost::get(mark_map, boost::target(*it, g->parent())))
                cnt += postorder(g, boost::target(*it, g->parent()), list + cnt, r);
        }
    }
    assert(boost::get(rank_map, g->local_to_global(v)) == r);
    list[cnt++] = v;
    return cnt;
}

// g: local
void Mincross::flat_reorder(graph_t * g)
{
    int i, j, r, pos, n_search, local_in_cnt, local_out_cnt, base_order;
    node_t v, t;
    node_t *left, *right;
    node_t *temprank = NULL;
    edge_t flat_e, e;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    if (!graph_propt[g].has_flat_edges)
        return;
    for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; ++r) {
        if (!graph_propt[g].rank[r].v.size()) continue;

        std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;

        base_order = boost::get(order_map, g->local_to_global(g_rank_r_v[0]));
        for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it)
            boost::put(mark_map, g->local_to_global(*it), false);

        temprank = new node_t[i+1];
        pos = 0;

        /* construct reverse topological sort order in temprank */
        std::vector<node_t>::iterator it = g_rank_r_v.begin();
        std::vector<node_t>::reverse_iterator rit = g_rank_r_v.rbegin();
        for (; it != g_rank_r_v.end(), rit != g_rank_r_v.rend(); ++it, ++rit) {
            v = (graph_propt[g].rankdir & 1) ? *it : *rit;

            local_in_cnt = local_out_cnt = 0;
            std::vector<edge_t> &flat_in_v = boost::get(flat_in_map, g->local_to_global(v));
            for (std::vector<edge_t>::iterator it_v = flat_in_v.begin(); it_v != flat_in_v.end(); ++it_v) {
            // assume constrain condition fullfilled
                if (constraining_flat_edge(g, v, *it_v)) local_in_cnt++;
            }
            std::vector<edge_t> &flat_out_v = boost::get(flat_out_map, g->local_to_global(v));
            for (std::vector<edge_t>::iterator it_v = flat_out_v.begin(); it_v != flat_out_v.end(); ++it_v) {
            // assume constrain condition fullfilled
                if (constraining_flat_edge(g, v, *it_v)) local_out_cnt++;
            }
            if ((local_in_cnt == 0) && (local_out_cnt == 0))
                temprank[pos++] = v;
            else {
                if (!boost::get(mark_map, g->local_to_global(v)) && local_in_cnt == 0) {
                    left = temprank + pos;
                    n_search = postorder(g, v, left, r);
                    pos += n_search;
                }
            }
        }

        if (pos) {
            if (!graph_propt[g].rankdir & 1) {
                left = temprank;
                right = temprank + pos - 1;
                while (left < right) {
                    t = *left;
                    *left = *right;
                    *right = t;
                    left++;
                    right--;
                }
            }

            i = 0;
            for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it, ++i) {
                v = (*it = temprank[i]);
                boost::put(order_map, g->local_to_global(v), i + base_order);
            }

            /* nonconstraint flat edges must be made LR */
            for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it) {
                std::vector<edge_t> flat_out_v = boost::get(flat_out_map, g->local_to_global(*it));
                if (flat_out_v.size()) {
                    for (std::vector<edge_t>::iterator it_v = flat_out_v.begin(); it_v != flat_out_v.end(); ++it_v) {
                        if (!(graph_propt[g].rankdir & 1) && boost::get(order_map, boost::target(*it_v, g->parent())) < 
                            boost::get(order_map, boost::source(*it_v, g->parent())) || (graph_propt[g].rankdir & 1) &&
                                boost::get(order_map, boost::target(*it_v, g->parent())) > 
                                    boost::get(order_map, boost::source(*it_v, g->parent()))) {
                            assert(!constraining_flat_edge(g, v, *it_v));

                            edge_t &e = *it_v;
                            if (boost::get(to_orig_map, *it_v) != edge_t() && 
                                boost::get(to_virt_map, boost::get(to_orig_map, *it_v)) == *it_v) 
                                    boost::put(to_virt_map, boost::get(to_orig_map, *it_v), edge_t());
                            it_v = flat_out_v.erase(it_v);
		                    --it_v;
                            std::vector<edge_t> &flat_in_e = boost::get(flat_in_map, boost::target(e, g->parent()));
                            for (std::vector<edge_t>::iterator it_in = flat_in_e.begin(); it_in != flat_in_e.end(); ++it_in)
                                if (*it_in == e) {
                                    flat_in_e.erase(it_in);
                                    break;
                                }
    
                            flat_rev(g, e);
                            boost::remove_edge(e, g->parent());
                        }
                    }
                }
            }
	    /* postprocess to restore intended order */
	    }
	    /* else do no harm! */
        //boost::ref_property_map<graph_t*, Agraphinfo_t>
        //        root_graph_propt(boost::get_property(*Root, graph_IDproperty));
        graph_propt[g].rank[r].valid = false;
    }
    if (temprank)
        free(temprank);
}

// g: local
void Mincross::reorder(graph_t * g, int r, int reverse, int hasfixed)
{
    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    int changed = 0, nelt;
    bool muststay;
    std::vector<node_t> &vlist = graph_propt[g].rank[r].v;    
    std::vector<node_t>::iterator lp, rp, ep = vlist.end();

    for (nelt = vlist.size() - 1; nelt >= 0; --nelt) {
        lp = vlist.begin();
        while (lp != ep) {
            /* find leftmost node that can be compared */
            while ((lp != ep) && (boost::get(mval_map, g->local_to_global(*lp)) < 0))
                ++lp;
            if (lp == ep)
                break;
            /* find the node that can be compared */
            muststay = false;
            rp = lp;
            ++rp;
            for (; rp != ep; ++rp) {
                // we don't have any cluster
                if (left2right(g, *lp, *rp)) {
                    muststay = true;
                    break;
                }
                if (boost::get(mval_map, g->local_to_global(*rp)) >= 0)
                    break;
                // we don't have any cluster
            }
            if (rp == ep)
                break;
            if (muststay == false) {
                register int p1 = boost::get(mval_map, g->local_to_global(*lp));
                register int p2 = boost::get(mval_map, g->local_to_global(*rp));
                if ((p1 > p2) || ((p1 == p2) && (reverse))) {
                    exchange(g, *lp, *rp);
                    changed++;
                }
            }
            lp = rp;
        }
        if ((hasfixed == false) && (reverse == false))
            --ep;
    }

    //boost::ref_property_map<graph_t*, Agraphinfo_t>
    //        root_graph_propt(boost::get_property(*Root, graph_IDproperty));

    if (changed) {
        graph_propt[g].rank[r].valid = false;
        if (r > 0)
            graph_propt[g].rank[r-1].valid = false;
    }
}

// g: local
void Mincross::mincross_step(graph_t * g, int pass)
{
    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    //boost::ref_property_map<graph_t*, Agraphinfo_t>
    //        root_graph_propt(boost::get_property(*Root, graph_IDproperty));

    int r, other, first, last, dir;
    bool hasfixed, reverse;

    if ((pass % 4) < 2)
        reverse = true;
    else
        reverse = false;
    if (pass % 2) {
        r = graph_propt[g].maxrank - 1;
        dir = -1;
    } /* up pass */
    else {
        r = 1;
        dir = 1;
    } /* down pass */

    if (pass % 2 == 0) {	/* down pass */
        first = graph_propt[g].minrank + 1;
        if (graph_propt[g].minrank > graph_propt[g].minrank)
            first--;
        last = graph_propt[g].maxrank;
        dir = 1;
    } else {			/* up pass */
        first = graph_propt[g].maxrank - 1;
        last = graph_propt[g].minrank;
        if (graph_propt[g].maxrank < graph_propt[g].maxrank)
            first++;
        dir = -1;
    }

    for (r = first; r != last + dir; r += dir) {
        other = r - dir;
        hasfixed = medians(g, r, other);
        reorder(g, r, reverse, hasfixed);
    }
    transpose(g, !reverse);
}

// g: local, e: global
void Mincross::flat_rev(graph_t * g, edge_t e)
{
    int j;
    edge_t rev;

    if (!boost::get(flat_out_map, boost::target(e, *Root)).size())
        rev = edge_t();
    else {
        std::vector<edge_t> &flat_out_head_e = boost::get(flat_out_map, boost::target(e, *Root));
        for (std::vector<edge_t>::iterator it = flat_out_head_e.begin(); it != flat_out_head_e.end(); ++it)
            if (boost::target((rev = *it), *Root) == boost::source(e, *Root))
                break;
    }

    if (rev != edge_t()) {
        fastgraph fg(g, propmap);
        fg.merge_oneway(e, rev);
        if (boost::get(to_virt_map, e) != edge_t())
            boost::put(to_virt_map, e, rev);

        if (boost::get(edge_type_map, rev) == FLATORDER
            && boost::get(to_orig_map, rev) == edge_t())
                boost::put(to_orig_map, rev, e);

        boost::get(other_map, boost::source(e, g->parent())).push_back(e);

    } else {
        fastgraph fg(g, propmap);
        rev = fg.new_virtual_edge(boost::target(e, g->parent()), boost::source(e, g->parent()), e);
        if (boost::get(edge_type_map, e) == FLATORDER)
            boost::put(edge_type_map, rev, FLATORDER);
        else
            boost::put(edge_type_map, rev, REVERSED);
        fg.flat_edge(&(g->parent()), rev);
    }
    
}

// g: local
void Mincross::flat_search(graph_t * g, node_t v)
{
    int i;
    bool hascl;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    adjmatrix_t *M = graph_propt[g].rank[boost::get(rank_map, g->local_to_global(v))].flat;

    boost::put(mark_map, g->local_to_global(v), true);
    boost::put(onstack_map, g->local_to_global(v), true);
    // our graph doesn't contain clusters
    if (boost::get(flat_out_map, g->local_to_global(v)).size()) {
        std::vector<edge_t> &flat_out_v = boost::get(flat_out_map, g->local_to_global(v));
        for (std::vector<edge_t>::iterator it = flat_out_v.begin(); it != flat_out_v.end(); ++it) {
            if (boost::get(weight_map, *it) == 0)
                continue;
            int flat_idx_head = boost::get(low_map, boost::target(*it, g->parent()));
            int flat_idx_tail = boost::get(low_map, boost::source(*it, g->parent()));
            if (boost::get(onstack_map, boost::target(*it, g->parent())) == true) {
                assert(flat_idx_head < M->nrows);
                assert(flat_idx_tail < M->ncols);
                M->data[flat_idx_head * M->ncols + flat_idx_tail] = 1;

                edge_t &e = *it;
                if (boost::get(to_orig_map, *it) != edge_t() && 
                    boost::get(to_virt_map, boost::get(to_orig_map, *it)) == *it) 
                        boost::put(to_virt_map, boost::get(to_orig_map, *it), edge_t());
                it = flat_out_v.erase(it);
		        --it;
                std::vector<edge_t> &flat_in_e = boost::get(flat_in_map, boost::target(e, g->parent()));
                for (std::vector<edge_t>::iterator it_in = flat_in_e.begin(); it_in != flat_in_e.end(); ++it_in)
                    if (*it_in == e) {
                        flat_in_e.erase(it_in);
                        break;
                    }

                if (boost::get(edge_type_map, e) == FLATORDER)
                    continue;
                flat_rev(g, e);
                boost::remove_edge(e, g->parent());                
            } else {
                assert(flat_idx_head < M->nrows);
                assert(flat_idx_tail < M->ncols);

                M->data[flat_idx_tail * M->ncols + flat_idx_head] = 1;
                if (!boost::get(mark_map, boost::target(*it, g->parent())))
                    flat_search(g, boost::target(*it, g->parent()));
            }
        }
    }
    boost::put(onstack_map, g->local_to_global(v), false);

}
// g: local
void Mincross::flat_breakcycles(graph_t * g)
{
    int i, r; 
    bool flat;
    node_t v;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));
    for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; ++r) {

        flat = false;

        std::vector<node_t> &rank_r_v = graph_propt[g].rank[r].v;
        i = 0;
        for (std::vector<node_t>::iterator it = rank_r_v.begin(); it != rank_r_v.end(); ++it, ++i) {
            boost::put(onstack_map, g->local_to_global(*it), false);
            boost::put(mark_map, g->local_to_global(*it), false);
            boost::put(low_map, g->local_to_global(*it), i);
            if (boost::get(flat_out_map, g->local_to_global(*it)).size() && !flat) {
                graph_propt[g].rank[r].flat = 
                    new_matrix(graph_propt[g].rank[r].v.size(), graph_propt[g].rank[r].v.size());
                flat = true;
            }
        }

        if (flat) {
            i = 0;
            for (std::vector<node_t>::iterator it = rank_r_v.begin(); it != rank_r_v.end(); ++it, ++i) {
                if (!boost::get(mark_map, g->local_to_global(*it)))  
                    flat_search(g, *it);
            }
        }
    }
}

adjmatrix_t* Mincross::new_matrix(int i, int j)
{
    adjmatrix_t *rv = new adjmatrix_t();
    rv->nrows = i;
    rv->ncols = j;
    rv->data = new bool[i*j];
    return rv;
}

struct less_than_edgeid {
    graph_t *Root;
    less_than_edgeid(graph_t *root) : Root(root) {}
    bool operator() (const edge_t &e0, const edge_t &e1)
    {
        return boost::get(boost::edge_index, *Root, e0) - boost::get(boost::edge_index, *Root, e1);
    }
};

int edgeidcmpf(int *id1, int *id2) {
    return (*id1) - (*id2);
}

// g: global
void Mincross::do_ordering_node (graph_t * g, node_t n, int outflag)
{
    int ne;
    node_t u, v;
    edge_t e, f, fe;

    edge_t *sortlist = TE_list;

    // assume there is no cluster
    ne = 0;
    if (outflag) {
        std::vector<edge_t> &out_n = boost::get(out_map, n);
        for (std::vector<edge_t>::iterator it = out_n.begin(); it != out_n.end(); ++it)
            // in our case, there is no between-cluster edges
            sortlist[ne++] = *it;
    } else {
        std::vector<edge_t> &in_n = boost::get(in_map, n);
        for (std::vector<edge_t>::iterator it = in_n.begin(); it != in_n.end(); ++it)
            // in our case, there is no between-cluster edges
            sortlist[ne++] = *it;
    }
    if (ne <= 1)
        return;
    /* write null terminator at end of list.
       requires +1 in TE_list alloccation */
    sortlist[ne] = edge_t();
    std::sort(sortlist, sortlist + ne, less_than_edgeid(Root));
    for (ne = 1; (f = sortlist[ne]) != edge_t(); ne++) {
        e = sortlist[ne - 1];
        if (outflag) {
            u = boost::target(e, *g);
            v = boost::target(f, *g);
        } else {
            u = boost::source(e, *g);
            v = boost::source(f, *g);
        }
        fastgraph fg(g, propmap);
        if (fg.find_flat_edge(g, u, v) != edge_t())
            return;
        fe = fg.new_virtual_edge(u, v, edge_t());
        boost::put(edge_type_map, fe, FLATORDER);
        fg.flat_edge(g, fe);
    }
}

void Mincross::do_ordering(graph_t * g, int outflag)
{
    /* Order all nodes in graph */
    node_it n, n_end;

    for (boost::tie(n, n_end) = boost::vertices(*g); n != n_end; ++n)
        do_ordering_node(g, *n, outflag);

}

/* merge connected components, create globally consistent rank lists */
// g: global
void Mincross::merge2(graph_t * g)
{
    int i, r;
    node_t v;

    /* merge the components and rank limits */

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    /* install complete ranks */
    graph_t::children_iterator gi, gi_end;
    for (boost::tie(gi, gi_end) = g->children(); gi != gi_end; ++gi) {
        boost::ref_property_map<graph_t*, Agraphinfo_t>
                local_graph_propt(boost::get_property(*gi, graph_IDproperty));

        for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; ++r) {
            std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
            std::vector<node_t> &local_rank_r_v = local_graph_propt[&(*gi)].rank[r].v;
            for (std::vector<node_t>::iterator it = local_rank_r_v.begin(); it != local_rank_r_v.end(); ++it)
                g_rank_r_v.push_back(gi->local_to_global(*it));
        }
    }

    for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; ++r) {
        std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
        i = 0;
        for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it, ++i) {
            if (*it == boost::graph_traits<graph_t>::null_vertex()) {
                if (Verbose)
                    fprintf(stderr,
	                    " rank %d has node %d < %d nodes\n",
	                        r, *it, g_rank_r_v.size()); // *it could be boost::get(vertex_index, *g, *it)

                break;
            }
            boost::put(order_map, *it, i);
        }
    }
}
// g: global
void Mincross::cleanup2(graph_t * g, int nc)
{
    int i, j, r, c;
    node_t v;
    edge_t e;

    if (TI_list) {
        delete[] TI_list;
        TI_list = NULL;
    }
    if (TE_list) {
        delete[] TE_list;
        TE_list = NULL;
    }

    /* remove node temporary edges for ordering nodes */
    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; ++r) {
        std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
        i = 0;
        for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it, ++i) {
	        boost::put(order_map, *it, i);
	        std::vector<edge_t> &flat_out_v = boost::get(flat_out_map, *it);
	        if (flat_out_v.size()) {
		        for (std::vector<edge_t>::iterator it_v = flat_out_v.begin(); it_v != flat_out_v.end(); ++it_v)
		            if (boost::get(edge_type_map, *it_v) == FLATORDER) {
		                fastgraph fg(g, propmap);
		                fg.delete_flat_edge(g, *it_v);
                        it_v--;
                    }
            }
        }
        adjmatrix_t* &matrix = graph_propt[g].rank[r].flat;
        if (matrix != NULL) {
            delete [] matrix->data;
            delete matrix;
        }
    }
    if (Verbose)
        fprintf(stderr, "mincross g: %d crossings, %.2f secs.\n",
            nc, elapsed_sec());
}

// g: global
void Mincross::dot_mincross(graph_t *g, bool doBalance) {

    int c, nc = 0;
    char *s;

    init_mincross(g);

    graph_t::children_iterator gi, gi_end;
    for (boost::tie(gi, gi_end) = g->children(); gi != gi_end; ++gi) {        // init_mccomp(&(*gi), c);
        nc += mincross(&(*gi), 0, 2, doBalance);
        dumpr(&(*gi));
    }
    merge2(g);
#ifdef DEBUG1
	dumpr();
#endif
    cleanup2(g, nc);
}
// g: local
void Mincross::dumpr(graph_t *g)
{
    int r;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    std::cout<<"check order in local graph:"<<g<<std::endl;
    for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; r++) {
        std::cout<<"Rank "<<r<<":"<<std::endl;
        std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
        for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it) {
            std::cout<<"Node "<<*it<<" has order:"<<boost::get(order_map, g->local_to_global(*it))<<std::endl;

        }
    }
}

void Mincross::dumpr()
{
    int r;
    graph_t *g = Root;

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    std::cout<<"check order:"<<std::endl;
    for (r = graph_propt[g].minrank; r <= graph_propt[g].maxrank; r++) {
        std::cout<<"Rank "<<r<<":"<<std::endl;
        std::vector<node_t> &g_rank_r_v = graph_propt[g].rank[r].v;
        for (std::vector<node_t>::iterator it = g_rank_r_v.begin(); it != g_rank_r_v.end(); ++it) {
            std::cout<<"Node "<<*it<<" has order:"<<boost::get(order_map, *it)<<std::endl;

        }
    }
}

