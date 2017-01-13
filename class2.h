#ifndef CLASS1_H
#define CLASS1_H

#include "types_bgl.h"
#include "propmap.h"

class Class2 {
    graph_t *G;
    PropertyMap propmap;
    orig_out_map_t *orig_out_map;
    out_map_t *out_map;
    rank_map_t *rank_map;
    to_virt_map_t *to_virt_map;
    weight_map_t *weight_map;
    xpenalty_map_t *xpenalty_map;
    count_map_t *count_map;
    edge_type_map_t *edge_type_map;
    /*/
    orig_in_map_t *orig_in_map;
    in_map_t *in_map;
    to_orig_map_t *to_orig_map;

//*/
    void make_chain(graph_t *g, node_t from, node_t to, edge_t orig);
    void merge_chain(graph_t *g, edge_t e, edge_t f, int flag);

public:
    Class2(graph_t *g, PropertyMap propmap) : G(g), propmap(propmap), orig_out_map(propmap.orig_out_map), out_map(propmap.out_map),
            rank_map(propmap.rank_map), to_virt_map(propmap.to_virt_map), weight_map(propmap.weight_map), xpenalty_map(propmap.xpenalty_map), 
            count_map(propmap.count_map), edge_type_map(propmap.edge_type_map) {}
    void class2(graph_t *g);

};

#endif // #ifndef CLASS1_H
