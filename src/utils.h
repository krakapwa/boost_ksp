/*
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributor(s): Laurent Lejeune (laurent.lejeune@artorg.unibe.ch).
 *
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "ksp.h"
#include "globals.h"

namespace utils {

    Vertex add_vertex(MyGraph & g, int id, std::string str);

    bool add_edge(MyGraph & g,
                     vertex_id_type n0,
                     vertex_id_type n1,
                     double w,
                     edge_id_type id,
                     std::string str_0,
                     std::string str_1,
                     int label);

    int get_max_id(const MyGraph & g);
    void print_edge(Edge e, const MyGraph & g);
    void print_dist_pred(std::vector<double> dist,
                        std::vector<Vertex> preds,
                        const MyGraph & g);
    std::tuple<VertexPath,bool> pred_to_path(std::vector<std::size_t> preds,
                                            const MyGraph & g, Vertex source,
                                            Vertex sink);
    void print_path(EdgeSet path, const MyGraph & g);
    void print_paths(EdgeSets paths, const MyGraph & g);
    void print_all(const MyGraph & g);

    void invert_edge(Edge e,
                     bool inv_label,
                     bool inv_algebraic_sign,
                     bool ignore_label_val,
                     MyGraph & g);

    void invert_edges(EdgeSet edge_path,
                      bool inv_label,
                      bool inv_algebraic_sign,
                      bool ignore_label_val,
                      MyGraph & g);


    void invert_edges(EdgeSets edge_sets,
                      bool inv_label,
                      bool inv_algebraic_sign,
                      bool ignore_label_val,
                      MyGraph & g);

    EdgeSets translate_edge_sets(EdgeSets P,
                                const MyGraph & g);

    Edge translate_edge(Edge e,
                        const MyGraph & g_p,
                        const MyGraph & g);

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

    void set_label_to_invalid_edges(EdgeSet e_in,
                                    MyGraph & g_in,
                                    MyGraph & g_out,
                                    int label,
                                    bool invert);
    void set_label_to_all(MyGraph & g, int label);
    void set_label(EdgeSet p, MyGraph & g, int label);

    EdgeSet remove_edge_from_set(Edge e,
                                 EdgeSet p,
                                 const MyGraph & g,
                                 bool inv_edge);

    EdgeSet remove_edge_from_set(Vertex v_in,
                                 Vertex v_out,
                                 EdgeSet p,
                                 const MyGraph & g,
                                 bool inv_edge);

    edge_id_type find_ind_edge_starting_with(EdgeSet p,
                                    Vertex v,
                                    const MyGraph & g);

    edge_id_type find_ind_edge_starting_with(EdgeSet p,
                                    vertex_id_type v_id,
                                    const MyGraph & g);

    edge_id_type find_ind_edge_ending_with(EdgeSet p,
                                  Vertex v,
                                  const MyGraph & g);

    double calc_cost(EdgeSets P, const MyGraph & g);
    Edge append_edge(EdgeSet p_app,
                     Vertex start,
                     const MyGraph & g);


    bp::list edgeSets_to_edges_list(EdgeSets P, const MyGraph & g);

    bp::list edgeSets_to_vertices_list(EdgeSets P, const MyGraph & g);

    EdgeSets augment2(EdgeSets P_l,
                     EdgeSet p_cut,
                     EdgeSet p_inter,
                     Vertex sink_vertex,
                     MyGraph & g);

    EdgeSets augment(EdgeSets P_l,
                     EdgeSet p_cut,
                     EdgeSet p_inter,
                     Vertex sink_vertex,
                     MyGraph & g);

    bool has_duplicate_vertex_ids(const MyGraph & g);
    bool has_duplicate_edge_ids(const MyGraph & g);
    bool has_duplicate_edge_ids(EdgeSet p, const MyGraph & g);
    int num_edges_with_label(const MyGraph & g, int label);
    EdgeSet out_edges_with_label(Vertex v, int label, const MyGraph & g);
    EdgeSet out_edges_with_neg_label(Vertex v, const MyGraph & g);
    EdgeSet out_edges_with_neg_label(Vertex v, EdgeSet& s, const MyGraph & g);
    Edge first_out_edge(Vertex v, EdgeSet& s,const MyGraph & g);

    Edge last_out_edge(Vertex v, EdgeSet& s,
                                           const MyGraph & g);
    EdgeSet flatten(EdgeSets & P, const MyGraph& g);

    EdgeSet out_edges_with_label(Vertex v, int label, EdgeSet & s,
                                 const MyGraph & g);
}

#endif
