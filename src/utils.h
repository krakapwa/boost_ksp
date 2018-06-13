#ifndef UTILS_H_
#define UTILS_H_

#include "pyksp.h"
#include "globals.h"

void print_edge(Edge e, const MyGraph & g);
void print_dist_pred(std::vector<double> dist,
                     std::vector<Vertex> preds,
                     const MyGraph & g);
VertexPath pred_to_path(std::vector<std::size_t> preds,
                        const MyGraph & g, Vertex source,
                        Vertex sink);
void print_path(VertexPath path, const MyGraph & g);
void print_path(EdgeSet path, const MyGraph & g);
void print_paths(EdgeSets paths, const MyGraph & g);
void print_all(const MyGraph & g);
EdgeSet vertpath_to_edgepath(VertexPath path, const MyGraph & g);

bool edge_is_in_set(Edge e, EdgeSet p, const MyGraph & g, bool inv_direction);

EdgeSet get_edges_from_label(const MyGraph & g, int label);
void set_label_to_edges(MyGraph & g, EdgeSet es, int label);

void set_label_to_all(MyGraph & g, int label);
EdgeSet remove_edge_from_set(Edge e,
                             EdgeSet p,
                             const MyGraph & g,
                             bool inv_edge);

int find_ind_edge_starting_with(EdgeSet p,
                                MyVertex v,
                                const MyGraph & g);

int find_ind_edge_ending_with(EdgeSet p,
                                MyVertex v,
                                const MyGraph & g);

std::tuple<EdgeSet, EdgeSet, EdgeSet, Edge>
append_inter(EdgeSet p,
             EdgeSet p_inter,
             EdgeSet p_cut,
             EdgeSet leftovers,
             MyVertex start,
             MyVertex sink,
             const MyGraph & g);

#endif
