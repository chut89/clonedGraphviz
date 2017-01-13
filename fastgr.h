#include "types_bgl.h"
#include "propmap.h"

class fastgraph {
    graph_t *G;
    
    orig_in_map_t *orig_in_map;
    orig_out_map_t *orig_out_map;
    in_map_t *in_map;
    out_map_t *out_map;
    to_virt_map_t *to_virt_map;
    to_orig_map_t *to_orig_map;
    flat_in_map_t *flat_in_map;
    flat_out_map_t *flat_out_map;
    node_type_map_t *node_type_map;

    minlen_map_t *minlen_map;
    weight_map_t *weight_map;
    count_map_t *count_map;
    xpenalty_map_t *xpenalty_map;
    edge_type_map_t *edge_type_map;

    edge_t ffe(node_t u, std::vector<edge_t> &ulist, node_t v, std::vector<edge_t> &vlist);
    void basic_merge(edge_t e, edge_t rep);

public:
    fastgraph(graph_t *g, PropertyMap propmap) : G(g), orig_in_map(propmap.orig_in_map), orig_out_map(propmap.orig_out_map),
        in_map(propmap.in_map), out_map(propmap.out_map), to_virt_map(propmap.to_virt_map), to_orig_map(propmap.to_orig_map),
        flat_in_map(propmap.flat_in_map), flat_out_map(propmap.flat_out_map), node_type_map(propmap.node_type_map),
        weight_map(propmap.weight_map), minlen_map(propmap.minlen_map), count_map(propmap.count_map), 
        xpenalty_map(propmap.xpenalty_map), edge_type_map(propmap.edge_type_map) {}

    node_t virtual_node();

    // fast graph
    edge_t new_virtual_edge(node_t u, node_t v, edge_t orig);
    edge_t find_fast_edge(node_t u, node_t v);
    void merge_oneway(edge_t e, edge_t rep);
    edge_t virtual_edge(node_t u, node_t v, edge_t orig);

    // flat edges
    void flat_edge(graph_t *g, edge_t e);
    edge_t find_flat_edge(graph_t *g, node_t u, node_t v);
    void delete_flat_edge(graph_t *g, edge_t e);

};
