#ifndef DOTINIT_H
#define DOTINIT_H

#include "types_bgl.h"
#include "propmap.h"

class dotinit {
    Agraphinfo_t graphinfo;
	rank_map_t rank_map;
	order_map_t order_map;
	orig_in_map_t orig_in_map;
	orig_out_map_t orig_out_map;
    in_map_t in_map;
    out_map_t out_map;
    to_virt_map_t to_virt_map;
    to_orig_map_t to_orig_map;
    node_type_map_t node_type_map;

	weight_map_t weight_map;
	minlen_map_t minlen_map;
    count_map_t count_map;
    xpenalty_map_t xpenalty_map;
    edge_type_map_t edge_type_map;

    PropertyMap propmap;

public:
    void dot_init_node(graph_t *g);
    void dot_init_edge(graph_t *g);
    void dot_init_node_edge(graph_t * g);
    void dot_layout(graph_t * g);
    void set_temp_orders(graph_t *g);
};

#endif

