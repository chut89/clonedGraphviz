#ifndef NS_H
#define NS_H

#include"types_bgl.h"
#include <setjmp.h>
#include "propmap.h"

extern void start_timer(void);
extern double elapsed_sec(void);

class NetworkSimplex {
	int init_graph(graph_t * g);
	void dfs_cutval(node_t v, edge_t par);
	int dfs_range(node_t v, edge_t par, int low);
	int x_val(edge_t e, node_t v, int dir);
	#ifdef DEBUG1
	void tchk();
	void check_cycles(graph_t * g);
	node_t checkdfs(graph_t *g, node_t n);
	void checktree(graph_t *g);
	int dump_node (graph_t *g, node_t n);
	void check_cutvalues();
	int check_ranks();

	#endif

    PropertyMap propmap;
	int length(edge_t edge);
	int slack(edge_t edge);
	//#define LENGTH(e)		(ND_rank(aghead(e)) - ND_rank(agtail(e)))
	//#define SLACK(e)		(LENGTH(e) - ED_minlen(e))
	//#define SEQ(a,b,c)		(((a) <= (b)) && ((b) <= (c)))
	//#define TREE_EDGE(e)	(boost::get(tree_index_map, e) >= 0)
	bool seq(int a, int b, int c);
	bool tree_edge(edge_t);

	jmp_buf jbuf;
	graph_t *G;
	int N_nodes, N_edges;
	int Minrank, Maxrank;
	int S_i;			/* search index for enter_edge */
	int Search_size;
	#define SEARCHSIZE 30
	#define	NORMAL 0	/* an original input node */
	bool Verbose;
	std::vector<node_t> Tree_node;
	std::vector<edge_t> Tree_edge;

	edge_t Enter;
	int Low, Lim, Slack;

	// property_maps
	mark_map_t mark_map;
	priority_map_t priority_map;
	rank_map_t rank_map;
	//order_map_t order_map;
	orig_in_map_t orig_in_map;
	orig_out_map_t orig_out_map;
    in_map_t in_map;
    out_map_t out_map;
	tree_in_map_t tree_in_map;
	tree_out_map_t tree_out_map;
	low_map_t low_map;
	lim_map_t lim_map;
	par_map_t par_map;
	node_type_map_t node_type_map;
	onstack_map_t onstack_map;

	cutvalue_map_t cutvalue_map;
	tree_index_map_t tree_index_map;
	minlen_map_t minlen_map;
	weight_map_t weight_map;
	//*/
	void add_tree_edge(edge_t e);
	void exchange_tree_edges(edge_t e, edge_t f);//*/
	void init_rank();
	node_t incident(edge_t e);
	edge_t leave_edge();
	void dfs_enter_outedge(node_t v);
	void dfs_enter_inedge(node_t v);
	edge_t enter_edge(edge_t e);
	bool treesearch(node_t v);
	int tight_tree();
	void init_cutvalues();
	int feasible_tree();
	node_t treeupdate(node_t v, node_t w, int cutvalue, bool dir);
	void rerank(node_t v, int delta);
	void update(edge_t e, edge_t f);
	void scan_and_normalize();
	void freeTreeList (graph_t *g);
	void LR_balance();
	void TB_balance();
	void graphSize (graph_t g, int* nn, int* ne);
	void x_cutval(edge_t f);
public:
	NetworkSimplex(graph_t *g, PropertyMap propmap);
    ~NetworkSimplex();
	int rank(graph_t *g, int balance, int maxiter);
	int rank2(graph_t *g, int balance, int maxiter, int search_size); 	

	void dump_graph (graph_t *g);

};

#endif // #define NS_H
