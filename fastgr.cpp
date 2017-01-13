#include "fastgr.h"
#include <iostream>

node_t fastgraph::virtual_node() {
    node_t n = boost::add_vertex(*G);
    boost::put(*node_type_map, n, VIRTUAL);
    //boost::put(boost::vertex_index, *G, n, boost::num_vertices(*G) - 1);
    return n;
}

edge_t fastgraph::find_fast_edge(node_t u, node_t v) {
    return ffe(u, boost::get(*out_map, u), v, boost::get(*in_map, v));
}

edge_t fastgraph::ffe(node_t u, std::vector<edge_t> &ulist, node_t v, std::vector<edge_t> &vlist) {
    if (ulist.size() < vlist.size()) 
        for (std::vector<edge_t>::iterator it = ulist.begin(); it != ulist.end(); ++it) {
            if (boost::target(*it, *G) == v)
                return *it;
        }
    else 
        for (std::vector<edge_t>::iterator it = vlist.begin(); it != vlist.end(); ++it) {
            if (boost::source(*it, *G) == u)
                return *it;
        }
    return edge_t();
}

edge_t fastgraph::virtual_edge(node_t u, node_t v, edge_t orig) {
    edge_t e = new_virtual_edge(u, v, orig);
    boost::get(*out_map, u).push_back(e);
    boost::get(*in_map, v).push_back(e);
    return e;
}

edge_t fastgraph::new_virtual_edge(node_t u, node_t v, edge_t orig) {
    edge_t e;
    bool inserted;
    boost::tie(e, inserted) = boost::add_edge(u, v, *G);

    if (!inserted) {
        std::cerr<<"Failed adding virtual edge"<<std::endl; 
        return edge_t();
    }
    boost::put(*edge_type_map, e, VIRTUAL);

    if (orig != edge_t()) {
        boost::put(*count_map, e, boost::get(*count_map, orig));
        boost::put(*xpenalty_map, e, boost::get(*xpenalty_map, orig));
	    boost::put(*weight_map, e, boost::get(*weight_map, orig));
        boost::put(*minlen_map, e, boost::get(*minlen_map, orig));

	    if (boost::get(*to_virt_map, orig) == edge_t())
	        boost::put(*to_virt_map, orig, e);
	    boost::put(*to_orig_map, e, orig);
    } else {
        boost::put(*count_map, e, 1);
        boost::put(*xpenalty_map, e, 1);
	    boost::put(*weight_map, e, 1);
        boost::put(*minlen_map, e, 1);
    }
    return e;
}

void fastgraph::basic_merge(edge_t e, edge_t rep) {
    if (boost::get(*minlen_map, rep) < boost::get(*minlen_map, e))
        boost::put(*minlen_map, rep, boost::get(*minlen_map, e));

    while (rep != edge_t()) { // in our case, there is only one level of virtual edge
        boost::put(*count_map, rep, boost::get(*count_map, rep) + boost::get(*count_map, e));
        boost::put(*xpenalty_map, rep, boost::get(*xpenalty_map, rep) + boost::get(*xpenalty_map, e));
        boost::put(*weight_map, rep, boost::get(*weight_map, rep) + boost::get(*weight_map, e));
        rep = boost::get(*to_virt_map, rep);
    }

}

void fastgraph::merge_oneway(edge_t e, edge_t rep) {
    if (rep == boost::get(*to_virt_map, e)) {
        std::cerr<< "merge_oneway glitch"<<std::endl;
        return;
    }
    assert(boost::get(*to_virt_map, e) == edge_t());

    boost::put(*to_virt_map, e, rep);
    basic_merge(e, rep);

}

edge_t fastgraph::find_flat_edge(graph_t *g, node_t u, node_t v) {
    return ffe(u, boost::get(*flat_out_map, u), v, boost::get(*flat_in_map, v));
}

void fastgraph::flat_edge(graph_t *g, edge_t e) {
    boost::get(*flat_out_map, boost::source(e, *g)).push_back(e);
    boost::get(*flat_out_map, boost::target(e, *g)).push_back(e);

    boost::ref_property_map<graph_t*, Agraphinfo_t>
            graph_propt(boost::get_property(*g, graph_IDproperty));

    graph_propt[g].has_flat_edges = true;
}

void fastgraph::delete_flat_edge(graph_t *g, edge_t e) {
    assert(e != edge_t());
    if (boost::get(*to_orig_map, e) != edge_t() && 
        boost::get(*to_virt_map, boost::get(*to_orig_map, e)) == e) 
            boost::put(*to_virt_map, boost::get(*to_orig_map, e), edge_t());

    std::vector<edge_t> flat_out_e = boost::get(*flat_out_map, boost::source(e, *g));
    for (std::vector<edge_t>::iterator it = flat_out_e.begin(); it != flat_out_e.end(); ++it)
        if (*it == e) {
            flat_out_e.erase(it);
            break;
        }
    std::vector<edge_t> flat_in_e = boost::get(*flat_in_map, boost::target(e, *g));
    for (std::vector<edge_t>::iterator it = flat_in_e.begin(); it != flat_in_e.end(); ++it)
        if (*it == e) {
            flat_in_e.erase(it);
            break;
        }

    boost::remove_edge(e, *g);
}

