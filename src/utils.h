#ifndef UTILS_H_
#define UTILS_H_

#include "pyksp.h"
#include "globals.h"

namespace utils {

    Vertex add_vertex(MyGraph & g, int id, std::string str);

    bool add_edge(MyGraph & g,
                     int n0,
                     int n1,
                     double w,
                     int id,
                     std::string str_0,
                     std::string str_1,
                     int label);
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

    void invert_edge(Edge e,
                     bool inv_label,
                     bool inv_algebraic_sign,
                     const MyGraph & g_in,
                     MyGraph & g_out);

    void invert_edges(EdgeSet edge_path,
                      bool inv_label,
                      bool inv_algebraic_sign,
                      const MyGraph & g_in,
                      MyGraph & g_out);


    void invert_edges(EdgeSets edge_sets,
                      bool inv_label,
                      bool inv_algebraic_sign,
                      const MyGraph & g_in,
                      MyGraph & g_out);
    EdgeSet translate_edge_set(EdgeSet p,
                               const MyGraph & g_p,
                               const MyGraph & g,
                               bool inv_mode);

    EdgeSet vertpath_to_edgepath(VertexPath path, const MyGraph & g);

    bool edge_is_in_set(Edge e,
                        EdgeSet p,
                        const MyGraph & g,
                        bool inv_direction);

    bool edge_is_in_set(Edge e,
                        EdgeSet p,
                        const MyGraph & g_e,
                        const MyGraph & g_p,
                        bool inv_direction);

    bool vertex_is_in_set(Vertex v,
                          EdgeSet p,
                          const MyGraph & g);

    EdgeSet get_edges_from_label(const MyGraph & g, int label);
    void set_label_to_edges(MyGraph & g, EdgeSet es, int label);

    void set_label_to_all(MyGraph & g, int label);
    void set_label_on_path(EdgeSet p,
                           int label,
                           MyGraph & g);

    EdgeSet remove_edge_from_set(Edge e,
                                 EdgeSet p,
                                 const MyGraph & g,
                                 bool inv_edge);

    EdgeSet remove_edge_from_set(Vertex v_in,
                                 Vertex v_out,
                                 EdgeSet p,
                                 const MyGraph & g,
                                 bool inv_edge);

    int find_ind_edge_starting_with(EdgeSet p,
                                    Vertex v,
                                    const MyGraph & g);

    int find_ind_edge_ending_with(EdgeSet p,
                                  Vertex v,
                                  const MyGraph & g);

    std::tuple<EdgeSet, EdgeSet, EdgeSet, Edge>
        append_inter(EdgeSet p,
                     EdgeSet p_inter,
                     EdgeSet p_cut,
                     EdgeSet leftovers,
                     Vertex start,
                     Vertex sink,
                     const MyGraph & g,
                     const MyGraph & g_c);
}

#endif
