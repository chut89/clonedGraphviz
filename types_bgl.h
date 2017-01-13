#ifndef TYPES_H
#define TYPES_H


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/incremental_components.hpp>
//#include <boost/graph/subgraph.hpp>
#include "subgraph.hpp"

enum graph_IDproperty_t
{
   graph_IDproperty
};

namespace boost
{
  BOOST_INSTALL_PROPERTY(graph, IDproperty);
}

struct Agnodeinfo_t;
struct Agedgeinfo_t;
struct Agraphinfo_t;

typedef boost::property<boost::vertex_index_t, std::size_t, Agnodeinfo_t> vertex_prop;
typedef boost::property<boost::edge_index_t, std::size_t, Agedgeinfo_t> edge_prop;
typedef boost::property<graph_IDproperty_t, Agraphinfo_t> graph_prop;
// listS as vertex container while being used with subgraph causes a lot of pains!
typedef boost::subgraph<boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS, vertex_prop, edge_prop, graph_prop> > graph_t;

//typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Agnodeinfo_t, Agedgeinfo_t> graph_t;//bgl
typedef boost::graph_traits<graph_t>::vertex_iterator node_it;//bgl
typedef boost::graph_traits<graph_t>::vertex_descriptor node_t;//bgl
typedef boost::graph_traits <graph_t>::vertices_size_type Size;//bgl
typedef boost::graph_traits<graph_t>::edge_iterator edge_it;//bgl
typedef boost::graph_traits<graph_t>::edge_descriptor edge_t;//bgl

typedef boost::graph_traits<graph_t>::in_edge_iterator in_edge_iterator;
typedef boost::graph_traits<graph_t>::out_edge_iterator out_edge_iterator;

typedef boost::subgraph<boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, vertex_prop, edge_prop> > undirected_graph_t;
//typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> undirected_graph_t;

//typedef boost::component_index<node_t> Components;

typedef struct pointf_s { double x, y; } pointf;

typedef struct Agnodeinfo_t {
	/*/
	Agrec_t hdr;
	shape_desc *shape;
	void *shape_info;//*/
	pointf coord;
	//double width, height;  /* inches */
	/*/
	boxf bb;
	double ht, lw, rw;
	textlabel_t *label;
	textlabel_t *xlabel;
	void *alg;
	char state;//*/
	//unsigned char gui_state; /* Node state for GUI ops */
	//boolean clustnode;
/*/
#ifndef DOT_ONLY
	unsigned char pinned;
	int id, heapindex, hops;
	double *pos, dist;
#endif//*/
//#ifndef NEATO_ONLY
	/*/
	unsigned char showboxes;
	boolean  has_port;
	node_t* rep;
	node_t *set;//*/

	/* fast graph */
	bool mark;
	char node_type;
	bool onstack;
	char ranktype, weight_class;
	node_t *next, *prev;
	std::vector<edge_t> orig_in, orig_out, in, out, flat_out, flat_in, other;//bgl
	//graph_t *clust;

	/* for union-find and collapsing nodes */
	int UF_size;
	node_t *UF_parent;//bgl
	node_t *inleaf, *outleaf;//bgl

	/* for placing nodes */
	int rank, order;	/* initially, order = 1 for ordered edges */
	double mval;
	std::vector<edge_t*> save_in, save_out;

	/* for network-simplex */
	std::vector<edge_t> tree_in, tree_out;//bgl
	edge_t par;//bgl
	int low, lim;
	int priority;

	double pad[1];
//#endif

} Agnodeinfo_t;

#define NORMAL 0
#define VIRTUAL 1
#define FLATORDER 2
#define REVERSED 3

typedef struct Agedgeinfo_t {
	/*/
	Agrec_t hdr;
	splines *spl;
	port tail_port, head_port;
	//textlabel_t *label, *head_label, *tail_label, *xlabel;//*/
	char edge_type;
	//char adjacent;          /* true for flat edge with adjacent nodes */
	//char label_ontop;
	//unsigned char gui_state; /* Edge state for GUI ops */
	edge_t to_orig;	/* for dot's shapes.c    */
	edge_t to_virt;
	//void *alg;
/*/
#ifndef DOT_ONLY
	double factor;
	double dist;
	//Ppolyline_t path;
#endif//*/
//#ifndef NEATO_ONLY
	//unsigned char showboxes;
	//boolean conc_opp_flag;
	short xpenalty;
	int weight;
	int cutvalue, tree_index;
	short count;
	unsigned short minlen;
//#endif
} Agedgeinfo_t;

typedef struct adjmatrix_t {
	int nrows, ncols;
	bool *data;
} adjmatrix_t;

typedef struct rank_t {
	//int n;			/* number of nodes in this rank  */
	std::vector<node_t> v;		/* ordered list of nodes in rank    */
	//int an;			/* globally allocated number of nodes   */
	//node_t **av;		/* allocated list of nodes in rank  */
	//double ht1, ht2;	/* height below/above centerline    */
	//double pht1, pht2;	/* as above, but only primitive nodes   */
	bool candidate;	/* for transpose () */
	bool valid;
	int cache_nc;		/* caches number of crossings */
	adjmatrix_t *flat;
} rank_t;

typedef struct Agraphinfo_t {
	//Agrec_t hdr;
	/* to generate code */
	//layout_t *drawing;
	//textlabel_t *label;	/* if the cluster has a title */
	//boxf bb;			/* bounding box */
	//pointf border[4];	/* sizes of margins for graph labels */
	//unsigned char gui_state; /* Graph state for GUI ops */
	//unsigned char has_labels;
	//boolean has_images;
	//unsigned char charset; /* input character set */
	int rankdir;
	double ht1, ht2;	/* below and above extremal ranks */
	unsigned short flags;
	//void *alg;
	//GVC_t *gvc;	/* context for "globals" over multiple graphs */
	//void (*cleanup) (graph_t * g);   /* function to deallocate layout-specific data */

//#ifndef DOT_ONLY
	/* to place nodes */
/*/
	node_t **neato_nlist;
	int move;
	double **dist, **spring, **sum_t, ***t;
	unsigned short ndim;
	unsigned short odim;
#endif//*/
//#ifndef NEATO_ONLY
	/* to have subgraphs */
//	int n_cluster;
//	graph_t **clust;	/* clusters are in clust[1..n_cluster] !!! */
	graph_t *dotroot;
//	node_it nlist;//bgl
	rank_t *rank;
//	graph_t *parent;        /* containing cluster (not parent subgraph) */
//	int level;		/* cluster nesting level (not node level!) */
	node_t	*minrep, *maxrep;	/* set leaders for min and max rank */

	/* fast graph node list */
//	std::vector<graph_t> comp;//bgl
	/* connected components */
	node_t *minset, *maxset;	/* set leaders */
	long n_nodes;
	/* includes virtual */
	short minrank, maxrank;

	/* various flags */
	bool has_flat_edges;
	bool has_sourcerank;
	bool has_sinkrank;
	//unsigned char	showboxes;
	//fontname_kind fontnames;		/* to override mangling in SVG */

	int nodesep, ranksep;
	//node_t *ln, *rn;	/* left, right nodes of bounding box */

	/* for clusters */
	/*/
	node_t *leader, **rankleader;
	boolean expanded;
	char installed;
	char set_type;
	char label_pos;
	boolean exact_ranksep;//*/
//#endif
	int searchsize;

} Agraphinfo_t;


typedef boost::property_map<graph_t, int Agnodeinfo_t::*>::type rank_map_t;
typedef boost::property_map<graph_t, int Agnodeinfo_t::*>::type order_map_t;
typedef boost::property_map<graph_t, bool Agnodeinfo_t::*>::type mark_map_t;
typedef boost::property_map<graph_t, int Agnodeinfo_t::*>::type priority_map_t;
typedef boost::property_map<graph_t, std::vector<edge_t> Agnodeinfo_t::*>::type orig_in_map_t;
typedef boost::property_map<graph_t, std::vector<edge_t> Agnodeinfo_t::*>::type orig_out_map_t;
typedef boost::property_map<graph_t, std::vector<edge_t> Agnodeinfo_t::*>::type in_map_t;
typedef boost::property_map<graph_t, std::vector<edge_t> Agnodeinfo_t::*>::type out_map_t;
typedef boost::property_map<graph_t, std::vector<edge_t> Agnodeinfo_t::*>::type tree_in_map_t;
typedef boost::property_map<graph_t, std::vector<edge_t> Agnodeinfo_t::*>::type tree_out_map_t;
typedef boost::property_map<graph_t, std::vector<edge_t> Agnodeinfo_t::*>::type flat_in_map_t;
typedef boost::property_map<graph_t, std::vector<edge_t> Agnodeinfo_t::*>::type flat_out_map_t;
typedef boost::property_map<graph_t, std::vector<edge_t> Agnodeinfo_t::*>::type other_map_t;
typedef boost::property_map<graph_t, int Agnodeinfo_t::*>::type low_map_t;
typedef boost::property_map<graph_t, int Agnodeinfo_t::*>::type lim_map_t;
typedef boost::property_map<graph_t, edge_t Agnodeinfo_t::*>::type par_map_t;
typedef boost::property_map<graph_t, char Agnodeinfo_t::*>::type node_type_map_t;
typedef boost::property_map<graph_t, bool Agnodeinfo_t::*>::type onstack_map_t;
typedef boost::property_map<graph_t, double Agnodeinfo_t::*>::type mval_map_t;
typedef boost::property_map<graph_t, pointf Agnodeinfo_t::*>::type coord_map_t;

typedef boost::property_map<graph_t, int Agedgeinfo_t::*>::type cutvalue_map_t;
typedef boost::property_map<graph_t, int Agedgeinfo_t::*>::type tree_index_map_t;
typedef boost::property_map<graph_t, unsigned short Agedgeinfo_t::*>::type minlen_map_t;
typedef boost::property_map<graph_t, int Agedgeinfo_t::*>::type weight_map_t;
typedef boost::property_map<graph_t, short Agedgeinfo_t::*>::type count_map_t;
typedef boost::property_map<graph_t, short Agedgeinfo_t::*>::type xpenalty_map_t;
typedef boost::property_map<graph_t, edge_t Agedgeinfo_t::*>::type to_virt_map_t;
typedef boost::property_map<graph_t, edge_t Agedgeinfo_t::*>::type to_orig_map_t;
typedef boost::property_map<graph_t, char Agedgeinfo_t::*>::type edge_type_map_t;

#endif // ifndef TYPES_H
