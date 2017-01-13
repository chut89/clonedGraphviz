#ifndef RANK_H
#define RANK_H

#include "types_bgl.h"
#include "ns.h"
#include "propmap.h"

class Rank {

	in_map_t *in_map;
	out_map_t *out_map;
    mark_map_t mark_map;
    rank_map_t *rank_map;
    to_virt_map_t *to_virt_map;
    to_orig_map_t *to_orig_map;
    edge_type_map_t *edge_type_map;
    void acyclic(graph_t *g);
	void rank1(graph_t *g);
    void expand_ranksets(graph_t * g);
    void cleanup1(graph_t * g);

	PropertyMap propmap;

public:
    Rank(PropertyMap propmap) : propmap(propmap), in_map(propmap.in_map), out_map(propmap.out_map),
        rank_map(propmap.rank_map), to_virt_map(propmap.to_virt_map), 
        to_orig_map(propmap.to_orig_map), edge_type_map(propmap.edge_type_map) {}
	void dot_rank(graph_t * g);
	
};

#endif
