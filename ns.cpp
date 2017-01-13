/* $Id$ $Revision$ */
/* vim:set shiftwidth=4 ts=8: */

/*************************************************************************
 * Copyright (c) 2011 AT&T Intellectual Property 
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors: See CVS logs. Details at http://www.graphviz.org/
 *************************************************************************/


/* 
 * Network Simplex Algorithm for Ranking Nodes of a DAG
 */

//#include "render.h"
//#include "memory.h"
//#include "cgraph.h"
//#include "globals.h"
//#include "const.h"
#include "ns.h"
#include <iostream>
#include <cstdio>
#include <queue>
#include <algorithm>

NetworkSimplex::NetworkSimplex(graph_t *g, PropertyMap propmap) : G(g), propmap(propmap) { 
//*/

    rank_map = boost::get(&Agnodeinfo_t::rank, *G);
    orig_in_map = boost::get(&Agnodeinfo_t::orig_in, *G);
    orig_out_map = boost::get(&Agnodeinfo_t::orig_out, *G);
    in_map = boost::get(&Agnodeinfo_t::in, *G);
    out_map = boost::get(&Agnodeinfo_t::out, *G);

    mark_map = boost::get(&Agnodeinfo_t::mark, *G);  
    priority_map = boost::get(&Agnodeinfo_t::priority, *G);
    tree_in_map = boost::get(&Agnodeinfo_t::tree_in, *G);
    tree_out_map = boost::get(&Agnodeinfo_t::tree_out, *G);
    par_map = boost::get(&Agnodeinfo_t::par, *G);
    low_map = boost::get(&Agnodeinfo_t::low, *G);
    lim_map = boost::get(&Agnodeinfo_t::lim, *G);
    node_type_map = boost::get(&Agnodeinfo_t::node_type, *G);
    onstack_map = boost::get(&Agnodeinfo_t::onstack, *G);

    node_it n, vi_end;
    for (boost::tie(n, vi_end) = boost::vertices(*G); n != vi_end; ++n) {
        boost::put(rank_map, *n, boost::get(*(propmap.rank_map), G->local_to_global(*n)));

        std::vector<edge_t> &orig_in_n = boost::get(*(propmap.orig_in_map), G->local_to_global(*n));
        std::vector<edge_t> orig_in_n_loc;
        for (std::vector<edge_t>::iterator it = orig_in_n.begin(); it != orig_in_n.end(); ++it)
            orig_in_n_loc.push_back(G->global_to_local(*it));
        boost::put(orig_in_map, *n, orig_in_n_loc);

        std::vector<edge_t> &orig_out_n = boost::get(*(propmap.orig_out_map), G->local_to_global(*n));
        std::vector<edge_t> orig_out_n_loc;
        for (std::vector<edge_t>::iterator it = orig_out_n.begin(); it != orig_out_n.end(); ++it)
            orig_out_n_loc.push_back(G->global_to_local(*it));
        boost::put(orig_out_map, *n, orig_out_n_loc);

        std::vector<edge_t> &in_n = boost::get(*(propmap.in_map), G->local_to_global(*n));
        std::vector<edge_t> in_n_loc;
        for (std::vector<edge_t>::iterator it = in_n.begin(); it != in_n.end(); ++it)
            in_n_loc.push_back(G->global_to_local(*it));
        boost::put(in_map, *n, in_n_loc);

        std::vector<edge_t> &out_n = boost::get(*(propmap.out_map), G->local_to_global(*n));
        std::vector<edge_t> out_n_loc;
        for (std::vector<edge_t>::iterator it = out_n.begin(); it != out_n.end(); ++it)
            out_n_loc.push_back(G->global_to_local(*it));
        boost::put(out_map, *n, out_n_loc);

        boost::put(mark_map, *n, false); 
        boost::put(priority_map, *n, 0);
        boost::put(par_map, *n, edge_t());
        boost::put(low_map, *n, 0);
        boost::put(lim_map, *n, 0);
        boost::put(node_type_map, *n, NORMAL);
        boost::put(onstack_map, *n, false);
    }

    minlen_map = boost::get(&Agedgeinfo_t::minlen, *G);
    weight_map = boost::get(&Agedgeinfo_t::weight, *G);

    cutvalue_map = boost::get(&Agedgeinfo_t::cutvalue, *G);
    tree_index_map = boost::get(&Agedgeinfo_t::tree_index, *G);

    edge_it e, e_end;
    for (boost::tie(e, e_end) = boost::edges(*G); e != e_end; ++e) {
        boost::put(minlen_map, *e, boost::get(*(propmap.minlen_map), G->local_to_global(*e)));
        boost::put(weight_map, *e, boost::get(*(propmap.weight_map), G->local_to_global(*e)));

        boost::put(cutvalue_map, *e, 0);
        boost::put(tree_index_map, *e, -1);
    }

    Enter = edge_t();

    Verbose = true;
//*/
    if (Verbose) checktree(g);
}

NetworkSimplex::~NetworkSimplex() {
    node_it n, vi_end;
    for (boost::tie(n, vi_end) = boost::vertices(*G); n != vi_end; ++n) {

        std::vector<edge_t> &orig_in_n = boost::get(orig_in_map, *n);
        std::vector<edge_t> orig_in_n_loc;
        for (std::vector<edge_t>::iterator it = orig_in_n.begin(); it != orig_in_n.end(); ++it)
            orig_in_n_loc.push_back(G->local_to_global(*it));
        boost::put(orig_in_map, *n, orig_in_n_loc);

        std::vector<edge_t> &orig_out_n = boost::get(orig_out_map, *n);
        std::vector<edge_t> orig_out_n_loc;
        for (std::vector<edge_t>::iterator it = orig_out_n.begin(); it != orig_out_n.end(); ++it)
            orig_out_n_loc.push_back(G->local_to_global(*it));
        boost::put(orig_out_map, *n, orig_out_n_loc);
    }    
}

void NetworkSimplex::add_tree_edge(edge_t e) {

    if (tree_edge(e)) {
	std::cerr<<"add_tree_edge: missing tree edge"<<std::endl;
	longjmp (jbuf, 1);
    }
    boost::put(tree_index_map, e, Tree_edge.size());
    Tree_edge.push_back(e);

    node_t tail = boost::source(e, *G);
    if (!boost::get(mark_map, tail)) 
	Tree_node.push_back(tail); 

    node_t head = boost::target(e, *G);
    if (!boost::get(mark_map, head)) 
	Tree_node.push_back(head); 

    boost::put(mark_map, tail, true);
    boost::get(tree_out_map, tail).push_back(e);
    
    boost::put(mark_map, head, true);
    boost::get(tree_in_map, head).push_back(e);    
    
}

void NetworkSimplex::exchange_tree_edges(edge_t e, edge_t f) {
    int i, j;
    node_t n;

    boost::put(tree_index_map, f, boost::get(tree_index_map, e));
    Tree_edge[boost::get(tree_index_map, e)] = f;
    boost::put(tree_index_map, e, -1);
    n = boost::source(e, *G);
    std::vector<edge_t> &tree_out_n = boost::get(tree_out_map, n);
    i = tree_out_n.size() - 1;
    for (j = 0; j <= i; j++)
	if (tree_out_n[j] == e) 
	    break;    
    tree_out_n[j] = tree_out_n[i];
    tree_out_n.pop_back();

    n = boost::target(e, *G);
    std::vector<edge_t> &tree_in_n = boost::get(tree_in_map, n);
    i = tree_in_n.size() - 1;
    for (j = 0; j <= i; j++)
	if (tree_in_n[j] == e)
	    break;
    tree_in_n[j] = tree_in_n[i];
    tree_in_n.pop_back();

    n = boost::source(f, *G);    
    boost::get(tree_out_map, n).push_back(f);

    n = boost::target(f, *G);
    boost::get(tree_in_map, n).push_back(f);

}//*/

void NetworkSimplex::init_rank() {
    int i, ctr;
    std::queue<node_t> Q;
    node_it v, vi_end;
    //in_edge_iterator ei, ei_end;
    //out_edge_iterator eo, eo_end;

    ctr = 0;

    for (boost::tie(v, vi_end) = boost::vertices(*G); v != vi_end; ++v) {
	if (boost::get(priority_map, *v) == 0)
	    Q.push(*v);
    }

    while (!Q.empty()) {
	node_t popped_node = Q.front();
	Q.pop();
	boost::put(rank_map, popped_node, 0);
	++ctr;
        std::vector<edge_t> in_popped_node = boost::get(in_map, popped_node);
        for (std::vector<edge_t>::iterator ei = in_popped_node.begin(); ei != in_popped_node.end(); ++ei)
	//for (boost::tie(ei, ei_end) = boost::in_edges(popped_node, *G); ei != ei_end; ++ei) 
	    boost::put(rank_map, popped_node, std::max(boost::get(rank_map, popped_node), 
		    boost::get(rank_map, boost::source(*ei, *G)) + boost::get(minlen_map, *ei)));
        
        std::vector<edge_t> out_popped_node = boost::get(out_map, popped_node);
        for (std::vector<edge_t>::iterator eo = out_popped_node.begin(); eo != out_popped_node.end(); ++eo) {
	//for (boost::tie(eo, eo_end) = boost::out_edges(popped_node, *G); eo != eo_end; ++eo) {
	    node_t target = boost::target(*eo, *G);
	    int pri_tmp = boost::get(priority_map, target);
	    --pri_tmp;
	    boost::put(priority_map, target, pri_tmp);
	    if (pri_tmp <= 0)
		Q.push(target);
	}
    }

    if (ctr != N_nodes) {
	std::cout<<"trouble in init_rank"<<std::endl;
    }
    
}

node_t NetworkSimplex::incident(edge_t e) {
    node_t tail = boost::source(e, *G);
    node_t head = boost::target(e, *G);
    if (boost::get(mark_map, tail)) {
	if (!boost::get(mark_map, head))
	    return tail;
    } else {
	if (boost::get(mark_map, head))
	    return head;
    }
    return boost::graph_traits<graph_t>::null_vertex();
}

edge_t NetworkSimplex::leave_edge() {
    edge_t f, rv = edge_t();
    int j, cnt = 0;

    j = S_i;
    while (S_i < Tree_edge.size()) {
        
        if (boost::get(cutvalue_map, f = Tree_edge[S_i]) < 0) {//cutvalues were wrongly calculated
	    if (rv != edge_t()) {
                if (boost::get(cutvalue_map, rv) > boost::get(cutvalue_map, f)) 
		    rv = f;
	    } else
		rv = Tree_edge[S_i];
	    if (++cnt >= Search_size)
		return rv;
	}
	S_i++;
    }
    if (j > 0) {
	S_i = 0;
	while (S_i < j) {
            if (boost::get(cutvalue_map, f = Tree_edge[S_i]) < 0) {
		if (rv != edge_t()) {
                    if (boost::get(cutvalue_map, rv) > boost::get(cutvalue_map, f))
			rv = f;
		} else
		    rv = Tree_edge[S_i];
		if (++cnt >= Search_size)
		    return rv;
	    }
	    S_i++;
	}
    } //std::cout<<"leave_edge: returned vertex="<<rv<<std::endl;
    return rv;
}//*/
//*/
void NetworkSimplex::dfs_enter_outedge(node_t v) {
    int i, slack_var;
    //out_edge_iterator eo, eo_end;

    std::vector<edge_t> out_v = boost::get(out_map, v);
    for (std::vector<edge_t>::iterator eo = out_v.begin(); eo != out_v.end(); ++eo) {
    //for (boost::tie(eo, eo_end) = boost::out_edges(v, *G); eo != eo_end; ++eo) {
        if (!tree_edge(*eo)) {
            if (!seq(Low, boost::get(lim_map, boost::target(*eo, *G)), Lim)) {
                slack_var = slack(*eo);
                if ((slack_var < Slack) || (Enter == edge_t())) { 
                    Enter = *eo;
                    Slack = slack_var;
                }		
            }
        } 
        else if (boost::get(lim_map, boost::target(*eo, *G)) < boost::get(lim_map, v))
            dfs_enter_outedge(boost::target(*eo, *G));	    
    }
    std::vector<edge_t> tree_in_tmp = boost::get(tree_in_map, v);
    for (std::vector<edge_t>::iterator it = tree_in_tmp.begin(); it != tree_in_tmp.end() && (Slack > 0); ++it) {
        node_t tail = boost::source(*it, *G);
        if (boost::get(lim_map, tail) < boost::get(lim_map, v))
            dfs_enter_outedge(tail);

    }

}
//*/
bool NetworkSimplex::seq(int a, int b, int c) {
    return (a <= b) && (b <= c);
}

bool NetworkSimplex::tree_edge(edge_t e) {
    return boost::get(tree_index_map, e) >= 0;
}
//*/
void NetworkSimplex::dfs_enter_inedge(node_t v) {
    int i, slack_var;
    //in_edge_iterator ei, ei_end;
    std::vector<edge_t> in_v = boost::get(in_map, v);
    for (std::vector<edge_t>::iterator ei = in_v.begin(); ei != in_v.end(); ++ei) {
    //for (boost::tie(ei, ei_end) = boost::in_edges(v, *G); ei != ei_end; ++ei) {
        if (!tree_edge(*ei)) {
            if (!seq(Low, boost::get(lim_map, boost::target(*ei, *G)), Lim)) {
                slack_var = slack(*ei);
                if ((slack_var < Slack) || (Enter == edge_t())) {
                    Enter = *ei;
                    Slack = slack_var;
                }		
            }
        } else if (boost::get(lim_map, boost::source(*ei, *G)) < boost::get(lim_map, v))
            dfs_enter_inedge(boost::source(*ei, *G));	    
    }
    std::vector<edge_t> tree_out_tmp = boost::get(tree_out_map, v);
    for (std::vector<edge_t>::iterator it = tree_out_tmp.begin(); it != tree_out_tmp.end() && (Slack > 0); ++it) {
        node_t head = boost::target(*it, *G);
        if (boost::get(lim_map, head) < boost::get(lim_map, v))
            dfs_enter_inedge(head);

    }
    
}//*/

edge_t NetworkSimplex::enter_edge(edge_t e) {
    node_t v;
    bool outsearch;

    /* v is the down node */
    node_t tail = boost::source(e, *G);
    node_t head = boost::target(e, *G);
    if (boost::get(lim_map, tail) < boost::get(lim_map, head)) {
        v = tail;
        outsearch = false;
    } else {
        v = head;
        outsearch = true;
    }
    
    Slack = INT_MAX;
    Low = boost::get(low_map, v);
    Lim = boost::get(lim_map, v);

    if (outsearch)
        dfs_enter_outedge(v);
    else
        dfs_enter_inedge(v);std::cout<<"enter_edge: returned vertex="<<v<<std::endl;
    return Enter;
}

bool NetworkSimplex::treesearch(node_t v) {
    int i;
    //in_edge_iterator ei, ei_end;
    //out_edge_iterator eo, eo_end;

    std::vector<edge_t> out_v = boost::get(out_map, v);
    for (std::vector<edge_t>::iterator eo = out_v.begin(); eo != out_v.end(); ++eo) {
        //for (boost::tie(eo, eo_end) = boost::out_edges(v, *G); eo != eo_end; ++eo) {
        if (!boost::get(mark_map, boost::target(*eo, *G)) && (slack(*eo) == 0)) {std::cout<<v<<" -> "<<boost::target(*eo, *G)<<std::endl;
            add_tree_edge(*eo); 
            if ((Tree_edge.size() == N_nodes - 1) || treesearch(boost::target(*eo, *G))) 
                return true;	    
        } 
    }

    std::vector<edge_t> in_v = boost::get(in_map, v);
    for (std::vector<edge_t>::iterator ei = in_v.begin(); ei != in_v.end(); ++ei) {
        //for (boost::tie(ei, ei_end) = boost::in_edges(v, *G); ei != ei_end; ++ei) {
        if (!boost::get(mark_map, boost::source(*ei, *G)) && (slack(*ei) == 0)) {std::cout<<boost::source(*ei, *G)<<" -> "<<v<<std::endl;
            add_tree_edge(*ei); 
            if ((Tree_edge.size() == N_nodes - 1) || treesearch(boost::source(*ei, *G)))
            return true;
        } 
    }
    return false;
}

int NetworkSimplex::length(edge_t e) {
    return boost::get(rank_map, boost::target(e, *G)) - boost::get(rank_map, boost::source(e, *G));
}
int NetworkSimplex::slack(edge_t e) {
    return length(e) - boost::get(minlen_map, e);
}

int NetworkSimplex::tight_tree() {
    int i;
    node_it n, vi_end;

    for (boost::tie(n, vi_end) = vertices(*G); n != vi_end; ++n) {
        boost::put(mark_map, *n, false);
        boost::get(tree_in_map, *n).clear();
        boost::get(tree_out_map, *n).clear();
    }    

    for (i = 0; i < Tree_edge.size(); i++)
    boost::put(tree_index_map, Tree_edge[i], -1);

    Tree_node.clear();
    Tree_edge.clear();

    for (boost::tie(n, vi_end) = vertices(*G); n != vi_end && Tree_edge.size() == 0; ++n) 
        treesearch(*n); 
    return Tree_node.size();
}

void NetworkSimplex::init_cutvalues() {
    node_it start = boost::vertices(*G).first;  
    dfs_range(*start, edge_t(), 1);
    dfs_cutval(*start, edge_t()); 
}

int NetworkSimplex::feasible_tree() {
    int i, delta;
    node_it n, n_end;
    //out_edge_iterator f, f_end;
    edge_t e;
    
    if (N_nodes <= 1)
    return 0;
	
    while (tight_tree() < N_nodes) {
    e = edge_t();
    for (boost::tie(n, n_end) = boost::vertices(*G); n != n_end; ++n) {
        std::vector<edge_t> out_n = boost::get(out_map, *n);
        for (std::vector<edge_t>::iterator f = out_n.begin(); f != out_n.end(); ++f) {
        //for (boost::tie(f, f_end) = boost::out_edges(*n, *G); f != f_end; ++f) {
            if (!tree_edge(*f) && (incident(*f) != boost::graph_traits<graph_t>::null_vertex()) && 
                            ((e == edge_t()) || (slack(*f) < slack(e)))) {
            //*/
                std::cout<<boost::source(*f, *G)<<" -> "<<boost::target(*f, *G)<<std::endl;
                std::cout<<tree_edge(*f)<<std::endl;
                std::cout<<(incident(*f) == boost::graph_traits<graph_t>::null_vertex())<<std::endl;
                std::cout<<(e == edge_t())<<std::endl;  
                std::cout<<slack(*f)<<" "<<std::endl;//*/
                    e = *f;
            }
        }
    }

    if (e != edge_t()) {
        delta = slack(e);
        if (delta) {
	    if (incident(e) == boost::target(e, *G))
	        delta = -delta;
        for (i = 0; i < Tree_node.size(); i++)
            boost::put(rank_map, Tree_node[i], boost::get(rank_map, Tree_node[i]) + delta);
            std::cout<<"shift rank map by delta="<<delta<<std::endl;
            check_ranks();
        }
    } else {
#ifdef DEBUG
	    std::cerr<<"not in tight tree:"<<std::endl;
    	    for (boost::tie(n, n_end) = boost::vertices(*G); n != n_end; ++n) {
		for (i = 0; i < Tree_node.size(); i++)
		    if (Tree_node[i] == *n)
			break;
		if (i >= Tree_node.size())
		    std::cerr<<"error with vertex "<<(*n)<<std::endl;
		
	    }
#endif
	    return 1;
	}
    }
    init_cutvalues(); check_cutvalues();
    return 0;
}

/* walk up from v to LCA(v,w), setting new cutvalues. */
node_t NetworkSimplex::treeupdate(node_t v, node_t w, int cutvalue, bool dir) {
    edge_t e;
    bool d;

    while (!seq(boost::get(low_map, v), boost::get(lim_map, w), boost::get(lim_map, v))) {
        e = boost::get(par_map, v);
        if (v == boost::source(e, *G))
            d = dir;
        else
            d = !dir;
        if (d)
            boost::put(cutvalue_map, e, boost::get(cutvalue_map, e) + cutvalue);
        else
            boost::put(cutvalue_map, e, boost::get(cutvalue_map, e) - cutvalue);
        if (boost::get(lim_map, boost::source(e, *G)) > boost::get(lim_map, boost::target(e, *G)))
            v = boost::source(e, *G);
        else
            v = boost::target(e, *G);

    }

    return v;
}

void NetworkSimplex::rerank(node_t v, int delta) {
    int i;
    edge_t e;

    boost::put(rank_map, v, boost::get(rank_map, v) - delta);
    std::vector<edge_t> tree_out_tmp = boost::get(tree_out_map, v);
    for (std::vector<edge_t>::iterator it = tree_out_tmp.begin(); it != tree_out_tmp.end(); ++it)
        if (*it != boost::get(par_map, v))
            rerank(boost::target(*it, *G), delta);

    std::vector<edge_t> tree_in_tmp = boost::get(tree_in_map, v);
    for (std::vector<edge_t>::iterator it = tree_in_tmp.begin(); it != tree_in_tmp.end(); ++it)
        if (*it != boost::get(par_map, v))
            rerank(boost::source(*it, *G), delta);

}

/* e is the tree edge that is leaving and f is the nontree edge that
 * is entering.  compute new cut values, ranks, and exchange e and f.
 */
//*/
void NetworkSimplex::update(edge_t e, edge_t f) {std::cout<<"update: e="<<e<<" and f="<<f<<std::endl;
    int cutvalue, delta;
    node_t lca;

    delta = slack(f);
    /* "for (v = in nodes in tail side of e) do ND_rank(v) -= delta;" */
    if (delta > 0) {
	int s;
        s = boost::get(tree_in_map, boost::source(e, *G)).size() + 
                boost::get(tree_out_map, boost::source(e, *G)).size();
	if (s == 1)
	    rerank(boost::source(e, *G), delta);
	else {
            s = boost::get(tree_in_map, boost::target(e, *G)).size() +
                boost::get(tree_out_map, boost::target(e, *G)).size();
	    if (s == 1)
		rerank(boost::target(e, *G), -delta);
	    else {
                if (boost::get(lim_map, boost::source(e, *G)) < boost::get(lim_map, boost::target(e, *G)))
		    rerank(boost::source(e, *G), delta);
		else
		    rerank(boost::target(e, *G), -delta);
	    }
	}
    }

    cutvalue = boost::get(cutvalue_map, e);
    lca = treeupdate(boost::source(f, *G), boost::target(f, *G), cutvalue, 1);
    //*/
    if (treeupdate(boost::target(f, *G), boost::source(f, *G), cutvalue, 0) != lca) {
	std::cerr<<"update: mismatched lca in treeupdates"<<std::endl;
	longjmp (jbuf, 1);
    }//*/
    boost::put(cutvalue_map, f, -cutvalue);
    boost::put(cutvalue_map, e, 0);    
    

    exchange_tree_edges(e, f); std::cout<<"lca="<<lca<<" par="<<boost::get(par_map, lca)<<" low="<<boost::get(low_map, lca)<<std::endl;
    dfs_range(lca, boost::get(par_map, lca), boost::get(low_map, lca));
}

void NetworkSimplex::scan_and_normalize() {
    node_it n, n_end;

    Minrank = INT_MAX;
    Maxrank = -INT_MAX;

    for (boost::tie(n, n_end) = boost::vertices(*G); n != n_end; ++n) {
	if (boost::get(node_type_map, *n) == NORMAL) {
	    int rank_n = boost::get(rank_map, *n);
	    Minrank = std::min(Minrank, rank_n);
	    Maxrank = std::max(Maxrank, rank_n);
	}
    }

    if (Minrank != 0) {
        for (boost::tie(n, n_end) = boost::vertices(*G); n != n_end; ++n)
            boost::put(rank_map, *n, boost::get(rank_map, *n) - Minrank);

	Maxrank -= Minrank;
	Minrank = 0;
    }
}

void NetworkSimplex::freeTreeList (graph_t *g) {
    node_it n, vi_end;
    for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) {
	boost::get(tree_in_map, *n).clear();
	boost::get(tree_out_map, *n).clear();
	boost::put(mark_map, *n, false);
    }
}

void NetworkSimplex::LR_balance() {
    int i, delta;
    edge_t e, f;

    for (i = 0; i < Tree_edge.size(); i++) {
	e = Tree_edge[i];

        if (boost::get(cutvalue_map, e) == 0) {
	    f = enter_edge(e);
	    if (f == edge_t())
		continue;
	    delta = slack(f);
	    if (delta <= 1)
		continue;

            if (boost::get(lim_map, boost::source(e, *G)) < boost::get(lim_map, boost::target(e, *G)))
                rerank(boost::source(e, *G), delta / 2);
	    else
		rerank(boost::target(e, *G), -delta / 2);
	}
    }
    freeTreeList (G);
}

void NetworkSimplex::TB_balance() {
    node_it n, n_end;
    //in_edge_iterator ei, ei_end;
    //out_edge_iterator eo, eo_end;
    int i, low, high, choice, *nrank;
    int inweight, outweight;

    scan_and_normalize();

    /* find nodes that are not tight and move to less populated ranks */
    nrank = new int[Maxrank + 1];
    for (i = 0; i <= Maxrank; i++)
	nrank[i] = 0;

    for (boost::tie(n, n_end) = boost::vertices(*G); n != n_end; ++n) 
        if (boost::get(node_type_map, *n) == NORMAL)
            nrank[boost::get(rank_map, *n)]++;

    for (boost::tie(n, n_end) = boost::vertices(*G); n != n_end; ++n) {
        if (boost::get(node_type_map, *n) != NORMAL)
            continue;
        
	inweight = outweight = 0;
	low = 0;
	high = Maxrank;

        std::vector<edge_t> in_n = boost::get(in_map, *n);
        for (std::vector<edge_t>::iterator ei = in_n.begin(); ei != in_n.end(); ++ei) {
        //for (boost::tie(ei, ei_end) = boost::in_edges(*n, *G); ei != ei_end; ++ei) {
	    inweight += boost::get(weight_map, *ei);
            low = std::max(low, boost::get(rank_map, boost::source(*ei, *G)) + boost::get(minlen_map, *ei));
        }

        std::vector<edge_t> out_n = boost::get(out_map, *n);
        for (std::vector<edge_t>::iterator eo = out_n.begin(); eo != out_n.end(); ++eo) {
        //for (boost::tie(eo, eo_end) = boost::out_edges(*n, *G); eo != eo_end; ++eo) {
	    outweight += boost::get(weight_map, *eo);
            high = std::min(high, boost::get(rank_map, boost::target(*eo, *G)) - boost::get(minlen_map, *eo));
        }        

	if (low < 0)//*/
	    low = 0;		/* vnodes can have ranks < 0 */
	if (inweight == outweight) {
	    choice = low;
	    for (i = low + 1; i <= high; i++)
		if (nrank[i] < nrank[choice])
		    choice = i;
            nrank[boost::get(rank_map, *n)]--;
	    nrank[choice]++;
            boost::put(rank_map, *n, choice);
	}
        boost::get(tree_in_map, *n).clear();
        boost::get(tree_out_map, *n).clear();
        boost::put(mark_map, *n, false);
    }

    delete [] nrank;
}

int NetworkSimplex::init_graph(graph_t *g) {

    bool feasible;
    node_it n, vi_end;
    //in_edge_iterator ei, ei_end;

    S_i = 0;
    N_edges = boost::num_edges(*g);
    N_nodes = boost::num_vertices(*g);
    for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) 
        boost::put(mark_map, *n, false);
    
    if (!Tree_node.size()) Tree_node = std::vector<node_t>();
    if (!Tree_edge.size()) Tree_edge = std::vector<edge_t>();

    feasible = true;
    for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) {
        boost::put(priority_map, *n, 0);
        std::vector<edge_t> in_n = boost::get(in_map, *n);
        for (std::vector<edge_t>::iterator ei = in_n.begin(); ei != in_n.end(); ++ei) {
            //for (boost::tie(ei, ei_end) = boost::in_edges(*n, *g); ei != ei_end; ++ei) {
            boost::put(priority_map, *n, boost::get(priority_map, *n) + 1);
            boost::put(cutvalue_map, *ei, 0);
            boost::put(tree_index_map, *ei, -1);
            if (feasible && (boost::get(rank_map, boost::target(*ei, *g)) - boost::get(rank_map, boost::source(*ei, *g)) < boost::get(minlen_map, *ei)))
                feasible = false;
        }

        boost::put(tree_in_map, *n, std::vector<edge_t>());
        boost::put(tree_out_map, *n, std::vector<edge_t>());

    }
    return feasible;
}

/* graphSize:
 * Compute no. of nodes and edges in the graph
 */
/*/
void NetworkSimplex::graphSize (graph_t g, int* nn, int* ne) {
    int i, nnodes, nedges;
    node_t *n;
    edge_t *e;
   
    nnodes = nedges = 0;
    for (n = GD_nlist(g); n; n = ND_next(n)) {
	nnodes++;
	for (i = 0; (e = ND_out(n).list[i]); i++) {
	    nedges++;
	}
    }
    *nn = nnodes;
    *ne = nedges;
}//*/

/* rank:
 * Apply network simplex to rank the nodes in a graph.
 * Uses ED_minlen as the internode constraint: if a->b with minlen=ml,
 * rank b - rank a >= ml.
 * Assumes the graph has the following additional structure:
 *   A list of all nodes, starting at GD_nlist, and linked using ND_next.
 *   Out and in edges lists stored in ND_out and ND_in, even if the node
 *  doesn't have any out or in edges.
 * The node rank values are stored in ND_rank.
 * Returns 0 if successful; returns 1 if `he graph was not connected;
 * returns 2 if something seriously wrong;
 */

int NetworkSimplex::rank2(graph_t *g, int balance, int maxiter, int search_size) {
    int iter = 0, feasible;
    char *ns = "network simplex: ";
    edge_t e, f;

#ifdef DEBUG
    check_cycles(g);
#endif
    //*/
    if (Verbose) {
	std::cerr<<boost::num_vertices(*g)<<" vertices and "<<boost::num_edges(*g)<<" edges and maxiter="<<maxiter<<" and balance="<<balance<<std::endl;
	start_timer();
    }//*/
    feasible = init_graph(g);std::cout<<feasible<<std::endl;
    if (!feasible)
	init_rank();
    if (maxiter <= 0) {
	freeTreeList (g);
	return 0;
    }

    if (search_size >= 0)
	Search_size = search_size;
    else
	Search_size = SEARCHSIZE;

    if (setjmp (jbuf)) {
	return 2;
    }

    if (feasible_tree()) {
	freeTreeList (g);
	return 1;
    } 
    while ((e = leave_edge()) != edge_t()) {
	f = enter_edge(e);
	update(e, f);
	iter++;
	if (Verbose && (iter % 100 == 0)) {
	    if (iter % 1000 == 100)
		fputs(ns, stderr);
	    fprintf(stderr, "%d ", iter);
	    if (iter % 1000 == 0)
		fputc('\n', stderr);
	}
	if (iter >= maxiter)
	    break;
    } //check_cutvalues();
    switch (balance) {
    case 1:
	TB_balance();
	break;
    case 2:
	LR_balance();
	break;
    default:
	scan_and_normalize();
	freeTreeList (G);
	break;
    }
    if (Verbose) {
	if (iter >= 100)
	    fputc('\n', stderr);
	fprintf(stderr, "%s%d nodes %d edges %d iter %.2f sec\n",
		ns, N_nodes, N_edges, iter, elapsed_sec());

	check_ranks();    
    }

    //node_it n, n_end;
    //for (boost::tie(n, n_end) = boost::vertices(*g); n != n_end; ++n)
    //    boost::put(*(propmap.rank_map), g->local_to_global(*n), boost::get(rank_map, *n));

    return 0;
}

int NetworkSimplex::rank(graph_t *g, int balance, int maxiter) {
    char *s;
    int search_size;

    search_size = SEARCHSIZE;

    return rank2 (g, balance, maxiter, search_size);
}

/* set cut value of f, assuming values of edges on one side were already set */
void NetworkSimplex::x_cutval(edge_t f) {
    node_t v;
    //in_edge_iterator ei, ei_end;
    //out_edge_iterator eo, eo_end;
    int i, sum, dir;

    /* set v to the node on the side of the edge already searched */
    if (boost::get(par_map, boost::source(f, *G)) == f) {
	v = boost::source(f, *G);
	dir = 1;
    } else {
	v = boost::target(f, *G);
	dir = -1;
    }

    sum = 0;

    std::vector<edge_t> out_v = boost::get(out_map, v);
    for (std::vector<edge_t>::iterator eo = out_v.begin(); eo != out_v.end(); ++eo) 
    //for (boost::tie(eo, eo_end) = boost::out_edges(v, *G); eo != eo_end; ++eo)
	sum += x_val(*eo, v, dir);
    std::vector<edge_t> in_v = boost::get(in_map, v);
    for (std::vector<edge_t>::iterator ei = in_v.begin(); ei != in_v.end(); ++ei) 
    //for (boost::tie(ei, ei_end) = boost::in_edges(v, *G); ei != ei_end; ++ei)
	sum += x_val(*ei, v, dir);    

    boost::put(cutvalue_map, f, sum);

}

int NetworkSimplex::x_val(edge_t e, node_t v, int dir) {
    node_t other;
    int d, rv, f;

    if (boost::source(e, *G) == v)
	other = boost::target(e, *G);
    else
	other = boost::source(e, *G);

    if (!seq(boost::get(low_map, v), boost::get(lim_map, other), boost::get(lim_map, v))) {
	f = 1;
	rv = boost::get(weight_map, e);
    } else {
	f = 0;
	if (tree_edge(e))
	    rv = boost::get(cutvalue_map, e);
	else
	    rv = 0;
	rv -= boost::get(weight_map, e);
    }

    if (dir > 0) {
	if (boost::target(e, *G) == v)
	    d = 1;
	else
	    d = -1;

    } else {
	if (boost::source(e, *G) == v)
	    d = 1;
	else
	    d = -1;
    }
    if (f)
	d = -d;
    if (d < 0)
	rv = -rv;
    return rv;
}

void NetworkSimplex::dfs_cutval(node_t v, edge_t par) {
    int i;

    std::vector<edge_t> tree_out_v = boost::get(tree_out_map, v);
    for (std::vector<edge_t>::iterator it = tree_out_v.begin(); it != tree_out_v.end(); ++it)
	if (*it != par)
	    dfs_cutval(boost::target(*it, *G), *it);

    std::vector<edge_t> tree_in_v = boost::get(tree_in_map, v);
    for (std::vector<edge_t>::iterator it = tree_in_v.begin(); it != tree_in_v.end(); ++it)
	if (*it != par)
	    dfs_cutval(boost::source(*it, *G), *it);

    if (par != edge_t())
	x_cutval(par);
}

int NetworkSimplex::dfs_range(node_t v, edge_t par, int low) {
    int i, lim;

    lim = low;

    boost::put(par_map, v, par);
    boost::put(low_map, v, low);

    std::vector<edge_t> tree_out_v = boost::get(tree_out_map, v);
    for (std::vector<edge_t>::iterator it = tree_out_v.begin(); it != tree_out_v.end(); ++it) {
        if (lim<20) std::cout<<"v="<<v<<" out="<<*it<<std::endl;
	if (*it != par)
	    lim = dfs_range(boost::target(*it, *G), *it, lim);  }
    std::vector<edge_t> tree_in_v = boost::get(tree_in_map, v);

    for (std::vector<edge_t>::iterator it = tree_in_v.begin(); it != tree_in_v.end(); ++it)
	if (*it != par)
	    lim = dfs_range(boost::source(*it, *G), *it, lim);  

    boost::put(lim_map, v, lim);
    if (lim<20) std::cout<<"end of dfs_range(v="<<v<<" par="<<par<<" lim="<<lim<<" low="<<low<<std::endl;
    return lim + 1;
}

#ifdef DEBUG1
void NetworkSimplex::tchk()
{
    int i, n_cnt, e_cnt;
    node_it n, n_end;
    edge_t *e;

    n_cnt = 0;
    e_cnt = 0;

    for (boost::tie(n, n_end) = boost::vertices(*G); n != n_end; ++n) {
        n_cnt++;
        std::vector<edge_t> tree_out_n = boost::get(tree_out_map, *n);
        for (std::vector<edge_t>::iterator it = tree_out_n.begin(); it != tree_out_n.end(); ++it) {
            e_cnt++;
            if (slack(*it) > 0)
            fprintf(stderr, "not a tight tree ");
        }
	
    }

    if ((n_cnt != Tree_node.size()) || (e_cnt != Tree_edge.size()))
        fprintf(stderr, "something missing\n");
}

void NetworkSimplex::check_cutvalues()
{
    node_it v, v_end;
    edge_t e;
    int i, save;

    for (boost::tie(v, v_end) = boost::vertices(*G); v != v_end; ++v) {
    std::vector<edge_t> tree_out_v = boost::get(tree_out_map, *v);
    for (std::vector<edge_t>::iterator it = tree_out_v.begin(); it != tree_out_v.end(); ++it) {
        save = boost::get(cutvalue_map, *it);
        x_cutval(*it);
        if (save != boost::get(cutvalue_map, *it))
        abort();std::cout<<"cutvalue:"<<save<<std::endl;
        }
    }
}

int NetworkSimplex::check_ranks()
{
    int cost = 0;
    node_it n, n_end;
    //out_edge_iterator eo, eo_end;

    for (boost::tie(n, n_end) = boost::vertices(*G); n != n_end; ++n) {
        std::vector<edge_t> out_n = boost::get(out_map, *n);
        for (std::vector<edge_t>::iterator eo = out_n.begin(); eo != out_n.end(); ++eo) {    
        //for (boost::tie(eo, eo_end) = boost::out_edges(*n, *G); eo != eo_end; ++eo) {
        cost += (boost::get(weight_map, *eo) * abs(length(*eo)));
        if (boost::get(rank_map, boost::target(*eo, *G)) - boost::get(rank_map, boost::source(*eo, *G)) - boost::get(minlen_map, *eo) < 0)
            abort();
        }
    }
    fprintf(stderr, "rank cost %d\n", cost);

    std::cout<<"rank list:"<<std::endl;
    for (boost::tie(n, n_end) = boost::vertices(*G); n != n_end; ++n) {
        std::cout<<boost::get(rank_map, *n)<<std::endl;
    }

    return cost;
}

void NetworkSimplex::checktree(graph_t *g)
{
    int i, n = 0, m = 0;
    node_it v, v_end;

    for (boost::tie(v, v_end) = boost::vertices(*g); v != v_end; ++v) {
        std::vector<edge_t> tree_out_v = boost::get(orig_out_map, *v);
        std::cout<<std::endl<<"node "<<(*v)<<std::endl<<"out edges: ";	
        for (std::vector<edge_t>::iterator it = tree_out_v.begin(); it != tree_out_v.end(); ++it) {
            std::cout<<boost::target(*it, *g)<<", ";
            ++n;
        }
        std::vector<edge_t> tree_in_v = boost::get(orig_in_map, *v);
        std::cout<<std::endl<<"in edges: ";
        for (std::vector<edge_t>::iterator it = tree_in_v.begin(); it != tree_in_v.end(); ++it) {
            std::cout<<boost::source(*it, *g)<<", ";
            ++m;
        }
    }

    for (boost::tie(v, v_end) = boost::vertices(*g); v != v_end; ++v) {
        std::vector<edge_t> tree_out_v = boost::get(out_map, *v);
        std::cout<<std::endl<<"node "<<(*v)<<std::endl<<"out edges: ";	
        for (std::vector<edge_t>::iterator it = tree_out_v.begin(); it != tree_out_v.end(); ++it) {
            std::cout<<boost::target(*it, *g)<<", ";
            ++n;
        }
        std::vector<edge_t> tree_in_v = boost::get(in_map, *v);
        std::cout<<std::endl<<"in edges: ";
        for (std::vector<edge_t>::iterator it = tree_in_v.begin(); it != tree_in_v.end(); ++it) {
            std::cout<<boost::source(*it, *g)<<", ";
            ++m;
        }
    }

    for (boost::tie(v, v_end) = boost::vertices(*g); v != v_end; ++v) {
        std::vector<edge_t> tree_out_v = boost::get(tree_out_map, *v);
        std::cout<<std::endl<<"node "<<(*v)<<std::endl<<"out edges: ";	
        for (std::vector<edge_t>::iterator it = tree_out_v.begin(); it != tree_out_v.end(); ++it) {
            std::cout<<boost::target(*it, *g)<<", ";
            ++n;
        }
        std::vector<edge_t> tree_in_v = boost::get(tree_in_map, *v);
        std::cout<<std::endl<<"in edges: ";
        for (std::vector<edge_t>::iterator it = tree_in_v.begin(); it != tree_in_v.end(); ++it) {
            std::cout<<boost::source(*it, *g)<<", ";
            ++m;
        }
    }

    fprintf(stderr, "%d %d %d\n", Tree_edge.size(), n, m);
}
/*/
void check_fast_node(node_it n)
{
    node_t *nptr;
    nptr = GD_nlist(agraphof(n));
    while (nptr && nptr != n)
	nptr = ND_next(nptr);
    assert(nptr != NULL);
}
//*/

int NetworkSimplex::dump_node (graph_t *g, node_t n)
{

    return static_cast<int>(boost::get(boost::vertex_index, *g, n));
}

void NetworkSimplex::dump_graph (graph_t *g)
{
    int i;
    edge_t *e;
    node_it n, n_end;
    //out_edge_iterator eo, eo_end;
    node_t w;
    FILE* fp = fopen ("ns.gv", "w");
    fprintf (fp, "digraph \"g\" {\n");

    for (boost::tie(n, n_end) = boost::vertices(*g); n != n_end; ++n) {
	fprintf (fp, "  \"%d\"\n", dump_node(g, *n));
        std::vector<edge_t> out_n = boost::get(out_map, *n);
        for (std::vector<edge_t>::iterator eo = out_n.begin(); eo != out_n.end(); ++eo) {
        //for (boost::tie(eo, eo_end) = boost::out_edges(*n, *g); eo != eo_end; ++eo) {
	    fprintf (fp, "  \"%d\"", dump_node(g, *n));
	    w = boost::target(*eo, *g);
	    fprintf (fp, " -> \"%d\"\n", dump_node(g, w));
	}
    }

    fprintf (fp, "}\n");
    fclose (fp);
}

node_t NetworkSimplex::checkdfs(graph_t *g, node_t n)
{
    //out_edge_iterator eo, eo_end;
    node_t w,x;
    int i;

    if (boost::get(mark_map, n))
	return boost::graph_traits<graph_t>::null_vertex();

    boost::put(mark_map, n, true);
    boost::put(onstack_map, n, true);

    std::vector<edge_t> out_n = boost::get(out_map, n);
    for (std::vector<edge_t>::iterator eo = out_n.begin(); eo != out_n.end(); ++eo) {
    //for (boost::tie(eo, eo_end) = boost::out_edges(n, *g); eo != eo_end; ++eo) {
	w = boost::target(*eo, *g);    

	if (boost::get(onstack_map, w)) {
	    dump_graph (g);
	    fprintf(stdout, "null edge: %d\n", (*eo == edge_t()));
	    fprintf(stderr, "cycle: last edge (%d) (%d)\n",
	       	 boost::get(boost::vertex_index, *g, n),
	       	 boost::get(boost::vertex_index, *g, w));
	    return w;
	}
	else {
	    if (!boost::get(mark_map, w)) {
		x = checkdfs(g, w);
		if (x != boost::graph_traits<graph_t>::null_vertex()) {
		    fprintf(stdout, "null edge: %d", (*eo == edge_t()));
		    fprintf(stderr,"unwind (%d)\n",
	         	 boost::get(boost::vertex_index, *g, n));
		    if (x != n) return x;
		    fprintf(stderr,"unwound to root\n");
		    fflush(stderr);
		    abort();
		    return boost::graph_traits<graph_t>::null_vertex();
		}
	    }
	}
    }
    boost::put(onstack_map, n, false);
    return boost::graph_traits<graph_t>::null_vertex();
}

void NetworkSimplex::check_cycles(graph_t *g)
{
    //*/
    node_it n, n_end;
    
    for (boost::tie(n, n_end) = boost::vertices(*g); n != n_end; ++n) {
	boost::put(mark_map, *n, false);
	boost::put(onstack_map, *n, false); 
    }

    for (boost::tie(n, n_end) = boost::vertices(*g); n != n_end; ++n) 
	checkdfs(g, *n);
    //*/
    dump_graph(g);
}
#endif				/* DEBUG */

