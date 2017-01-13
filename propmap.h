#ifndef PROPMAP_H
#define PROPMAP_H

class PropertyMap {

public:
	rank_map_t *rank_map;
	order_map_t *order_map;
	orig_in_map_t *orig_in_map;
	orig_out_map_t *orig_out_map;
    in_map_t *in_map;
    out_map_t *out_map;
    to_virt_map_t *to_virt_map;
    to_orig_map_t *to_orig_map;
    flat_in_map_t *flat_in_map;
    flat_out_map_t *flat_out_map;
    node_type_map_t *node_type_map;

	weight_map_t *weight_map;
	minlen_map_t *minlen_map;
    count_map_t *count_map;
    xpenalty_map_t *xpenalty_map;
    edge_type_map_t *edge_type_map;

    std::vector<int> components;
};

#endif 
