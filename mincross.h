#ifndef MINCROSS_H
#define MINCROSS_H

#include "types_bgl.h"
#include "propmap.h"
#include <queue>

extern void start_timer(void);
extern double elapsed_sec(void);

class Mincross {
    graph_t *Root;
    //int GlobalMinRank, GlobalMaxRank;

    // property maps
    PropertyMap propmap;

	rank_map_t rank_map;
	order_map_t order_map;
	orig_in_map_t orig_in_map;
	orig_out_map_t orig_out_map;
    in_map_t in_map;
    out_map_t out_map;
    flat_in_map_t flat_in_map;
    flat_out_map_t flat_out_map;
    other_map_t other_map;

	mark_map_t mark_map;
	low_map_t low_map;
	node_type_map_t node_type_map;
	onstack_map_t onstack_map;
	mval_map_t mval_map;
	coord_map_t coord_map;

	weight_map_t weight_map;
    xpenalty_map_t xpenalty_map;
    edge_type_map_t edge_type_map;
    to_orig_map_t to_orig_map;
    to_virt_map_t to_virt_map;

    #define		MC_SCALE	256
    int *TI_list;
    edge_t *TE_list;
    int MaxIter;
    double Convergence;
    bool Verbose;
    int MinQuit;
    bool ReMincross;

    typedef int(* qsort_cmpf)(const void *, const void *);
    void mincross_options(graph_t * g);
    void allocate_ranks(graph_t * g);
    void init_mincross(graph_t *g); 
    int ncross(graph_t * g);
    int rcross(graph_t * g, int r);
    bool left2right(graph_t * g, node_t v, node_t w);
    int in_cross(node_t v, node_t w);
    int out_cross(node_t v, node_t w);
    void exchange(graph_t *g, node_t v, node_t w);
    int transpose_step(graph_t * g, int r, bool reverse);
    void transpose(graph_t * g, bool reverse);
    void enqueue_neighbors(graph_t *g, std::queue<node_t> &q, node_t n0, int pass);
    void install_in_rank(graph_t * g, node_t n);
    void build_ranks(graph_t * g, int pass);
    bool flat_mval(graph_t *g, node_t n);
    bool medians(graph_t * g, int r0, int r1);
    bool is_a_normal_node_of(graph_t *g, node_t v);
    bool is_a_vnode_of_an_edge_of(graph_t *g, node_t v);
    bool inside_cluster(graph_t * g, node_t v);
    bool constraining_flat_edge(graph_t *g, node_t v, edge_t e);
    int postorder(graph_t * g, node_t v, node_t *list, int r);
    void flat_reorder(graph_t * g);
    void reorder(graph_t * g, int r, int reverse, int hasfixed);
    void balanceNodes(graph_t * g, int r, node_t v, node_t w);
    int balance(graph_t * g);
    void mincross_step(graph_t * g, int pass);
    int mincross(graph_t * g, int startpass, int endpass, int doBalance);
    void restore_best(graph_t * g);
    void save_best(graph_t * g);
    void flat_rev(graph_t * g, edge_t e);
    void flat_search(graph_t * g, node_t v);
    void flat_breakcycles(graph_t * g);
    adjmatrix_t *new_matrix(int i, int j);
    void do_ordering_node (graph_t * g, node_t n, int outflag);
    void do_ordering(graph_t * g, int outflag);
    void do_ordering_for_nodes(graph_t * g);
    void merge_components(graph_t * g);
    void merge2(graph_t *g);
    void cleanup2(graph_t * g, int nc);
#ifdef DEBUG1
    void dumpr();
    void dumpr(graph_t *g);
#endif

public:
    Mincross(graph_t *g, PropertyMap propmap);
    void dot_mincross(graph_t *g, bool doBalance);

};

#endif

