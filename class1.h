#ifndef CLASS1_H
#define CLASS1_H

#include "types_bgl.h"
#include "propmap.h"

class Class1 {
    graph_t *G;
    PropertyMap propmap;
    orig_in_map_t *orig_in_map;
    orig_out_map_t *orig_out_map;
    to_virt_map_t *to_virt_map;
    edge_type_map_t *edge_type_map;

public:
    Class1(graph_t *g, PropertyMap propmap) : G(g), propmap(propmap), orig_in_map(propmap.orig_in_map), orig_out_map(propmap.orig_out_map),
        to_virt_map(propmap.to_virt_map), edge_type_map(propmap.edge_type_map) {}
    void class1(graph_t *g);
};

#endif // #ifndef CLASS1_H
