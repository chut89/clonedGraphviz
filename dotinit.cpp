#include "dotinit.h"
#include "rank.h"
#include "mincross.h"
#include <iostream>

void dotinit::dot_init_node(graph_t *g) {
    node_it n, vi_end;
    for (boost::tie(n, vi_end) = boost::vertices(*g); n != vi_end; ++n) {
        boost::put(&Agnodeinfo_t::rank, *g, *n, 0);
        boost::put(&Agnodeinfo_t::order, *g, *n, 0);
        boost::put(&Agnodeinfo_t::orig_in, *g, *n, std::vector<edge_t>());
        boost::put(&Agnodeinfo_t::orig_out, *g, *n, std::vector<edge_t>());
        boost::put(&Agnodeinfo_t::in, *g, *n, std::vector<edge_t>());
        boost::put(&Agnodeinfo_t::out, *g, *n, std::vector<edge_t>());
        boost::put(&Agnodeinfo_t::node_type, *g, *n, NORMAL);

    }

    rank_map = boost::get(&Agnodeinfo_t::rank, *g);
    order_map = boost::get(&Agnodeinfo_t::order, *g);
    orig_in_map = boost::get(&Agnodeinfo_t::orig_in, *g);
    orig_out_map = boost::get(&Agnodeinfo_t::orig_out, *g);
    in_map = boost::get(&Agnodeinfo_t::in, *g);
    out_map = boost::get(&Agnodeinfo_t::out, *g);
    node_type_map = boost::get(&Agnodeinfo_t::node_type, *g);
    to_virt_map = boost::get(&Agedgeinfo_t::to_virt, *g);
    to_orig_map = boost::get(&Agedgeinfo_t::to_orig, *g);


    propmap.rank_map = &rank_map;
    propmap.order_map = &order_map;    
    propmap.orig_in_map = &orig_in_map;
    propmap.orig_out_map = &orig_out_map;
    propmap.in_map = &in_map;
    propmap.out_map = &out_map;
    propmap.node_type_map = &node_type_map;
    propmap.to_virt_map = &to_virt_map;
    propmap.to_orig_map = &to_orig_map;

}

void dotinit::dot_init_edge(graph_t *g) {

    edge_it e, e_end;
    for (boost::tie(e, e_end) = boost::edges(*g); e != e_end; ++e) {
        boost::put(&Agedgeinfo_t::to_virt, *g, *e, edge_t());
        boost::put(&Agedgeinfo_t::to_orig, *g, *e, edge_t());
        boost::put(&Agedgeinfo_t::minlen, *g, *e, 1);
        boost::put(&Agedgeinfo_t::weight, *g, *e, 1);
        boost::put(&Agedgeinfo_t::xpenalty, *g, *e, 1);
        boost::put(&Agedgeinfo_t::count, *g, *e, 1);
        boost::put(&Agedgeinfo_t::edge_type, *g, *e, NORMAL);
    }

    minlen_map = boost::get(&Agedgeinfo_t::minlen, *g);
    weight_map = boost::get(&Agedgeinfo_t::weight, *g);
    count_map = boost::get(&Agedgeinfo_t::count, *g);
    xpenalty_map = boost::get(&Agedgeinfo_t::xpenalty, *g);
    edge_type_map = boost::get(&Agedgeinfo_t::edge_type, *g);

    propmap.minlen_map = &minlen_map;
    propmap.weight_map = &weight_map;
    propmap.count_map = &count_map;
    propmap.xpenalty_map = &xpenalty_map;
    propmap.edge_type_map = &edge_type_map;

}

void dotinit::dot_init_node_edge(graph_t * g) {
    dot_init_node(g);
    dot_init_edge(g);

}

void dotinit::dot_layout(graph_t * g) {
    dot_init_node_edge(g);
    propmap.components = std::vector<int>(boost::num_vertices(*g));
    Rank r(propmap);
    r.dot_rank(g);
    Mincross mc(g, propmap);
    mc.dot_mincross(g, 0);
    //set_temp_orders(g);

}

void dotinit::set_temp_orders(graph_t *g) {
  // create arbirtrary order of nodes having same rank
  int max_rank = 0;
  node_it n, n_end;
  for (boost::tie(n, n_end) = boost::vertices(*g); n != n_end; ++n) { std::cout<<boost::get(rank_map, *n)<<std::endl;
    if (boost::get(rank_map, *n) > max_rank)
      max_rank = boost::get(rank_map, *n);
  } 
  int order[max_rank];
  for (int i = 0; i <= max_rank; ++i)
    order[i] = 0;
  for (boost::tie(n, n_end) = boost::vertices(*g); n != n_end; ++n) {
    boost::put(order_map, *n, order[boost::get(rank_map, *n)]++);
    //std::cout<<order[boost::get(rank_map, *n)]<<std::endl;
  }
}

