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

#include "utils.h"

namespace utils {
VertexDesc add_vertex(MyGraph &g, int id, std::string str) {
  // This has to pass through all vertices before adding (overhead).

  std::pair<VertexIter, VertexIter> vp;
  for (vp = vertices(g); vp.first != vp.second; ++vp.first)
    if (g[*(vp.first)].id == id) {
      // vertex was already there
      return *(vp.first);
    }

  // add new vertex
  VertexDesc res = boost::add_vertex(g);
  g[res].name = str;
  g[res].id = id;

  return res;
}

bool add_edge(MyGraph &g, vertex_id_type n0, vertex_id_type n1, double w,
              edge_id_type id = -1, std::string str_0 = "",
              std::string str_1 = "", int label = 1) {

  // Add two vertices
  VertexDesc v0 = add_vertex(g, n0, str_0);
  VertexDesc v1 = add_vertex(g, n1, str_1);

  std::pair<EdgeDesc, bool> e;
  e = boost::edge(v0, v1, g);

  if (!e.second) {
    e = boost::add_edge(v0, v1, g);

    g[e.first].weight = w;
    g[e.first].label = label;
    g[e.first].id = id;
    g[e.first].id_vertex_in = n0;
    g[e.first].id_vertex_out = n1;
    return true;
  } else
    return false;
}

std::tuple<VertexPath, bool> pred_to_path(std::vector<std::size_t> preds,
                                          const MyGraph &g, VertexDesc source,
                                          VertexDesc sink) {
  // Converts from predecessor map to vector of edge descriptors

  auto v_index = get(boost::vertex_index, g);

  VertexPath path;
  VertexDesc current = sink;
  VertexDesc last;
  bool ok = true;

  path.push_back(current);

  while (g[current].id != g[source].id) {

    last = current;
    current = preds[v_index[current]];

    if (g[current].id == g[last].id) {
      ok = false;
      break;
    }

    path.push_back(current);
  }

  // path.push_back(source);

  return std::make_tuple(path, ok);
}

void print_edge(const Edge& e, const MyGraph &g) {
  BOOST_LOG_TRIVIAL(trace) << "(" << g[e.in].name << ","
                           << g[e.out].name << ")"
                           << "/"
                           << "(" << g[e.in].id << ","
                           << g[e.out].id << ") "
                           << "label: " << g[e.to_edge_desc(g).first].label
                           << " weight: " << g[e.to_edge_desc(g).first].weight;
}

void print_path(const EdgeSet & path, const MyGraph &g) {

  for (unsigned int i = 0; i < path.size(); ++i) {
    print_edge(path[i], g);
  }
}

void print_paths(const EdgeSets & paths, const MyGraph &g) {

  for (unsigned int i = 0; i < paths.size(); ++i) {
    BOOST_LOG_TRIVIAL(trace) << "path #" << i;
    print_path(paths[i], g);
  }
}

EdgeSet get_edges_from_label(const MyGraph &g, int label) {
  EdgeIter ei, ei_end;
  EdgeSet out(g);

  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
    if (g[*ei].label == label)
      out += Edge(*ei, g);
  }
  return out;
}

void set_label(EdgeSet p, MyGraph &g, int label) {

  EdgeSet::iterator it;
  for (it = p.begin(); it != p.end(); ++it) {
    g[(*it).to_edge_desc(g).first].label = label;
  }
}

void set_label(EdgeDesc e, MyGraph &g, int label) { g[e].label = label; }

void set_label_to_all(MyGraph &g, int label) {
  // Assigns a label to all edges

  EdgeIter ei, ei_end;

  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
    g[*ei].label = label;
  }
}

// Assigns a label to edges that are invalid
void set_label_to_invalid_edges(EdgeSet p, MyGraph &g_out,
                                int label, bool invert) {
  // Assigns a label to invalid edges

  VertexDesc u;
  VertexDesc v;

  for (unsigned int i = 0; i < p.size(); ++i) {
    if (invert) {
      u = p[i].out;
      v = p[i].in;
    } else {
      u = p[i].in;
      v = p[i].out;
    }

    if (!edge(u, v, g_out).second)
      BOOST_LOG_TRIVIAL(trace)
          << "set_label_to_invalid_edges: Edge doesnt exist!!!";

    g_out[edge(u, v, g_out).first].label = label;
  }
}

void set_label_to_edges(MyGraph &g, EdgeSet p, int label) {
  // Assigns a label to edges

  for (unsigned int i = 0; i < p.size(); ++i) {
    g[p[i].to_edge_desc(g).first].label = label;
  }
}

void print_all(const MyGraph &g) {
  // Print the whole graph with node/edges ids, weights and labels

  BOOST_LOG_TRIVIAL(trace) << "------ graph: " << g[graph_bundle].name
                           << " --------";
  BOOST_LOG_TRIVIAL(trace) << "------ edges "
                           << " (n_edges: " << num_edges(g) << ") --------";
  EdgeIter ei, ei_end;
  MyGraph::vertex_descriptor v_in;
  MyGraph::vertex_descriptor v_out;

  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
    v_in = source(*ei, g);
    v_out = target(*ei, g);
    BOOST_LOG_TRIVIAL(trace)
        << "desc : " << *ei << " edge: (" << g[v_in].id << "," << g[v_out].id
        << ") / (" << g[v_in].name << "," << g[v_out].name
        << ") , id: " << g[*ei].id << ", label: " << g[*ei].label
        << ", weight: " << g[*ei].weight;
  }

  BOOST_LOG_TRIVIAL(trace) << "------ vertices"
                           << " (n_vertices: " << num_vertices(g)
                           << ") --------";
  std::pair<VertexIter, VertexIter> vp;
  for (vp = vertices(g); vp.first != vp.second; ++vp.first)
    BOOST_LOG_TRIVIAL(trace)
        << "desc: " << *vp.first << " id: " << g[*vp.first].id
        << ", name: " << g[*vp.first].name;
  BOOST_LOG_TRIVIAL(trace) << "--------------";
}

// Converts from vertex to edge
EdgeSet vertpath_to_edgepath(VertexPath v_path, const MyGraph &g) {
  EdgeSet p(g);
  // Convert set of vertices to set of edges

  VertexPath::reverse_iterator it;
  for (it = v_path.rbegin(); it != v_path.rend() - 1; ++it) {
    Edge e = Edge(*it, *std::next(it));
    p += e;
  }
  return p;
}

bool edge_is_in_set(EdgeDesc e, EdgeSet p, const MyGraph &g,
                    bool inv_direction = false) {
  // Is edge e in edge set p?

  VertexDesc curr_v_in;
  VertexDesc curr_v_out;

  VertexDesc v_in;
  VertexDesc v_out;

  if (inv_direction) {
    v_in = target(e, g);
    v_out = source(e, g);
  } else {
    v_in = source(e, g);
    v_out = target(e, g);
  }

  for (unsigned int i = 0; i < p.size(); ++i) {
    curr_v_in = p[i].in;
    curr_v_out = p[i].out;
    if ((g[curr_v_in].id == g[v_in].id) && (g[curr_v_out].id == g[v_out].id))
      return true;
  }

  return false;
}

bool edge_is_in_set(EdgeDesc e, EdgeSet p, const MyGraph &g_e, const MyGraph &g_p,
                    bool inv_direction = false) {
  // Is edge e in edge set p?

  VertexDesc v_in_p;
  VertexDesc v_out_p;

  VertexDesc v_in_e;
  VertexDesc v_out_e;

  if (inv_direction) {
    v_in_e = target(e, g_e);
    v_out_e = source(e, g_e);
  } else {
    v_in_e = source(e, g_e);
    v_out_e = target(e, g_e);
  }

  for (unsigned int i = 0; i < p.size(); ++i) {
    v_in_p = p[i].in;
    v_out_p = p[i].out;
    if ((g_p[v_in_p].id == g_e[v_in_e].id) &&
        (g_p[v_out_p].id == g_e[v_out_e].id))
      return true;
  }

  return false;
}

edge_id_type find_ind_edge_ending_with(EdgeSet p, VertexDesc v, const MyGraph &g) {
  // Search p for an edge ending with vertex v

  for (unsigned int i = 0; i < p.size(); ++i) {
    if (g[p[i].out].id == g[v].id)
      return i;
  }

  return -1;
}

edge_id_type find_ind_edge_starting_with(EdgeSet p, vertex_id_type v_id,
                                         const MyGraph &g) {
  // Search p for an edge starting with vertex v with id v_id

  BOOST_LOG_TRIVIAL(trace) << "find_ind_edge_starting_with (int): " << v_id;
  for (unsigned int i = 0; i < p.size(); ++i) {
    if (g[p[i].in].id == v_id) {
      return i;
    }
  }
  return -1;
}

edge_id_type find_ind_edge_starting_with(EdgeSet p, VertexDesc v,
                                         const MyGraph &g) {
  // Search p for an edge starting with vertex v

  for (unsigned int i = 0; i < p.size(); ++i) {
    if (g[p[i].in].id == g[v].id) {
      return i;
    }
  }

  return -1;
}

Edge append_edge(EdgeSet p_app, VertexDesc start, const MyGraph &g) {
  // This function is called when interlacing occurs.
  // One edge of p_app is returned

  // Find index of start in p_app
  edge_id_type ind_start = find_ind_edge_starting_with(p_app, g[start].id, g);

  return p_app[ind_start];
}

void invert_edge(const Edge & e, bool inv_label, bool inv_algebraic_sign,
                 bool ignore_label_val, MyGraph &g) {
  // invert edge e. Can invert sign of labels and/or
  // algebraic sign of weight

  VertexDesc u = e.in;
  VertexDesc v = e.out;

  // print_edge(e, g);

  if (g[e.to_edge_desc(g).first].label > 0 || ignore_label_val) {
    double weight = g[e.to_edge_desc(g).first].weight;
    int label = g[e.to_edge_desc(g).first].label;
    edge_id_type id = g[e.to_edge_desc(g).first].id;

    if (inv_algebraic_sign)
      weight = -weight;

    if (inv_label)
      label = -label;

    boost::remove_edge(u, v, g);

    add_edge(g, g[v].id, g[u].id, weight, id, g[v].name, g[u].name, label);
  }
}

void invert_edges(const EdgeSets & edge_sets,
                  bool inv_label, bool inv_algebraic_sign,
                  bool ignore_label_val, MyGraph &g) {
  // invert edges on edge_sets. Can invert sign of labels and/or
  // algebraic sign of weight

  for (unsigned int i = 0; i < edge_sets.size(); ++i) {
    invert_edges(edge_sets[i], inv_label, inv_algebraic_sign, ignore_label_val,
                 g);
  }
}

void invert_edges(const EdgeSet & edge_path, bool inv_label, bool inv_algebraic_sign,
                  bool ignore_label_val, MyGraph &g) {
  // invert path edge_path. Can invert sign of labels and/or
  // algebraic sign of weight

  for (unsigned int i = 0; i < edge_path.size(); ++i) {
    invert_edge(edge_path[i], inv_label, inv_algebraic_sign, ignore_label_val,
                g);
  }
}

double calc_cost(const EdgeSets & P, const MyGraph &g) {
  // Returns the total cost of edge sets

  double cost = 0;

  for (unsigned int i = 0; i < P.size(); ++i) {
    for (unsigned int j = 0; j < P[i].size(); ++j) {
      // BOOST_LOG_TRIVIAL(debug) << "calc_cost loop";
      cost += g[P[i][j].to_edge_desc(g).first].weight;
    }
  }

  return cost;
}

bp::list edgeSets_to_edges_list(EdgeSets P, const MyGraph &g) {
  // Convert to python list with edges ids

  bp::list array;
  for (unsigned int i = 0; i < P.size(); i++) {
    boost::python::list temp;
    for (unsigned int j = 0; j < P[i].size(); j++) {
      temp.append(g[P[i][j].to_edge_desc(g).first].id);
    }
    array.append(temp);
  }
  BOOST_LOG_TRIVIAL(debug) << "finished converting to list";

  return array;
}

bp::list edgeSets_to_vertices_list(EdgeSets P, const MyGraph &g) {
  // Convert to python list with edges ids

  bp::list array;
  VertexDesc u;
  VertexDesc v;
  for (unsigned int i = 0; i < P.size(); i++) {
    boost::python::list temp;
    for (unsigned int j = 0; j < P[i].size(); j++) {
      // temp.append(g[P[i][j]].id);
      u = P[i][j].in;
      v = P[i][j].out;
      // std::cout << "u: " << g[u].id << ", v: " << g[v].id << std::endl;
      temp.append(bp::make_tuple<int, int>(g[u].id, g[v].id));
    }
    array.append(temp);
  }

  return array;
}

EdgeSets augment(EdgeSets P_l, EdgeSet p_cut, EdgeSet p_inter,
                 VertexDesc sink_vertex, MyGraph &g) {
  /*
    P_l: Solution of past iteration
    p_cut: Edges of p_inter on inverted edges (are removed from candidates
    p_inter: Interlacing path

    Approach:
    - Start by augmenting p_inter using edges of P_l
    - Edges added to augmented path are labeled 0
    - Iterate (reverse order) through paths of P_l
    - If no edges have been taken previously, add path as is
    - Else start from first edge (never already occupied!)
        - Add edges until sink. Can choose to (1) take from remains
            of p_inter, or (2) take from other paths of P_l
   */

  // initialize P_l_plus_1 with empty sets
  EdgeSets P_l_plus_1(g);
  EdgeSet P_l_plus_1_flat(g);

  // This stores edges of last run to pick from
  EdgeSet P_l_clean = flatten(P_l, g) - p_cut;
  P_l_clean.sort_descend_labels();
  p_inter -= p_cut;

  BOOST_LOG_TRIVIAL(trace) << "Building p_ from p_inter and last solution";
  EdgeSet p_(g);

  print_path(p_inter, g);

  Edge curr_edge = p_inter[0];
  VertexDesc curr_vertex;
  p_ += curr_edge;
  print_edge(curr_edge, g);
  while (curr_edge.out != sink_vertex) {
    curr_vertex = curr_edge.out;
    BOOST_LOG_TRIVIAL(trace) << "curr_vertex: " << curr_vertex;
    BOOST_LOG_TRIVIAL(trace) << p_inter.has_out_vertex(curr_vertex);
    if (p_inter.has_out_vertex(curr_vertex)) {
      BOOST_LOG_TRIVIAL(trace) << "append inter";
      curr_edge = append_edge(p_inter, curr_vertex, g);
      // remove all edges that have been interlaced
      BOOST_LOG_TRIVIAL(trace) << "remove interlaced edges";

      // clean p_inter
      p_inter -= curr_edge;

    } else {
      BOOST_LOG_TRIVIAL(trace) << "pick from last solution set";
      curr_edge = Edge(last_out_edge(curr_vertex, P_l_clean, g), g);
      P_l_clean -= curr_edge;
    }
    p_ += curr_edge;
    curr_edge = p_.back();
    print_edge(curr_edge, g);
  }

  P_l_clean.sort_descend_labels();

  set_label(p_, g, 0);

  print_all(g);
  int label_p;
  EdgeSet edges_self(g);
  EdgeSets::reverse_iterator p;
  for (p = P_l.rbegin(); p != P_l.rend(); ++p) {
    if (!((*p).has_label(0))) { // ok (not cut by p_inter)
      BOOST_LOG_TRIVIAL(trace) << "inserting complete path:";
      print_path(*p, g);
      set_label(*p, g, 0);
      P_l_plus_1.insert(*p);
    } else {
      // we are tracking this label
      label_p = g[(*p)[0].to_edge_desc(g).first].label; 
      Edge curr_edge = (*p)[0];
      EdgeSet p_new(g);
      p_new += curr_edge;
      P_l_clean -= curr_edge;
      while (curr_edge.out != sink_vertex) {
        BOOST_LOG_TRIVIAL(trace)
            << "are P_l_clean labels sorted: " << P_l_clean.are_label_sorted();
        BOOST_LOG_TRIVIAL(trace) << "label: " << label_p;
        curr_vertex = curr_edge.out;
        edges_self = out_edges_with_label(curr_vertex, label_p, P_l_clean, g);
        // append inter or branch to other path
        if (edges_self.size() == 0) {
          if (p_inter.has_out_vertex(curr_vertex)) {
            BOOST_LOG_TRIVIAL(trace) << "append inter";
            curr_edge = append_edge(p_inter, curr_vertex, g);
            // remove all edges that have been interlaced
            BOOST_LOG_TRIVIAL(trace) << "remove interlaced edges";

            // clean p_inter
            p_new += curr_edge;
            p_inter -= curr_edge;
            P_l_clean -= curr_edge;
          } else {
            BOOST_LOG_TRIVIAL(trace) << "branch to another path";

            curr_edge = Edge(last_out_edge(curr_vertex, P_l_clean, g), g);
            p_new += curr_edge;
            P_l_clean -= curr_edge;
            BOOST_LOG_TRIVIAL(trace) << "end branch to another path";
          }
          curr_edge = p_new.back();
          // set_label(curr_edge, g, label_p);
        } else { // follow myself
          curr_edge = edges_self[0];
          p_new += curr_edge;
          P_l_clean -= curr_edge;
        }

        // print_edge(curr_edge, g);
      }
      BOOST_LOG_TRIVIAL(trace) << "lala";
      p_new.convert_to_graph(g);
      set_label(p_new, g, 0);
      BOOST_LOG_TRIVIAL(trace) << "lolo";
      P_l_plus_1.insert(p_new);
    }
  }

  BOOST_LOG_TRIVIAL(trace) << "done augmenting l paths";
  P_l_plus_1 += p_;

  BOOST_LOG_TRIVIAL(trace) << "settings labels to all";
  for (unsigned int i = 0; i < P_l_plus_1.size(); ++i)
    set_label(P_l_plus_1[i], g, -(i + 1));

  return P_l_plus_1;
}

bool has_duplicate_vertex_ids(const MyGraph &g) {

  std::vector<int> ids;
  VertexIter vi, vi_end;

  for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
    ids.push_back(g[*vi].id);
  }

  auto it = std::unique(ids.begin(), ids.end());
  return !(it == ids.end());
}

bool has_duplicate_edge_ids(EdgeSet p, const MyGraph &g) {

  std::vector<int> ids;
  EdgeSet::iterator ei;

  for (ei = p.begin(); ei != p.end(); ++ei) {
    ids.push_back(g[(*ei).to_edge_desc(g).first].id);
  }

  auto it = std::unique(ids.begin(), ids.end());
  return !(it == ids.end());
}

bool has_duplicate_edge_ids(const MyGraph &g) {

  std::vector<int> ids;
  EdgeIter ei, ei_end;

  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
    ids.push_back(g[*ei].id);
  }

  auto it = std::unique(ids.begin(), ids.end());
  return !(it == ids.end());
}

int num_edges_with_label(const MyGraph &g, int label) {

  std::vector<int> ids;
  EdgeIter ei, ei_end;

  for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
    if (g[*ei].label == label)
      ids.push_back(g[*ei].id);
  }

  return ids.size();
}

EdgeSet out_edges_with_label(VertexDesc v, int label, const MyGraph &g) {
  // return out_edges on graph g at vertex v with label label.

  EdgeSet out(g);
  MyGraph::out_edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
    if (g[*ei].label == label)
      out += Edge((*ei), g);
  }

  return out;
}

EdgeSet out_edges_with_label(VertexDesc v, int label, EdgeSet &p,
                             const MyGraph &g) {
  // return out_edges on graph g at vertex v with label label.

  EdgeSet out(g);
  EdgeSet::iterator it;
  for (it = p.begin(); it != p.end(); ++it) {
    if ((g[(*it).to_edge_desc(g).first].label == label) && ((*it).in == v))
      out += *it;
  }

  return out;
}

EdgeSet out_edges_with_neg_label(VertexDesc v, const MyGraph &g) {
  // return out_edges on graph g at vertex v with negative label.

  EdgeSet out(g);
  MyGraph::out_edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
    if (g[*ei].label < 0) {
      out += Edge(*ei, g);
    }
  }

  return out;
}

EdgeDesc first_out_edge(VertexDesc v, EdgeSet &p, const MyGraph &g) {
  // return edges on graph g at vertex v. Returns first instance.
  EdgeSet::iterator it;
  for (it = p.begin(); it != p.end(); ++it) {
    // g[(*it).in]
    if (g[(*it).in].id == g[v].id) {
      return (*it).to_edge_desc(g).first;
    }
  }

  BOOST_LOG_TRIVIAL(trace) << "first_out_edge: found no edge with in vertex: "
                           << g[v].id;
  return (*it).to_edge_desc(g).first;
}

EdgeDesc last_out_edge(VertexDesc v, EdgeSet &p, const MyGraph &g) {
  // return edges on graph g at vertex v. Returns last instance.
  EdgeSet::reverse_iterator it;
  for (it = p.rbegin(); it != p.rend(); ++it) {
    if (g[(*it).in].id == g[v].id) {
      BOOST_LOG_TRIVIAL(trace)
        << "last_out_edge: return label: " <<
        g[(*it).to_edge_desc(g).first].label;
      return (*it).to_edge_desc(g).first;
    }
  }
  BOOST_LOG_TRIVIAL(trace) << "last_out_edge: found no edge with in vertex: "
                           << g[v].id;

  return (*it).to_edge_desc(g).first;
}

EdgeSet out_edges_with_neg_label(VertexDesc v, EdgeSet &p, const MyGraph &g) {
  // return edges of set s on graph g at vertex v with negative label.

  EdgeSet out(g);
  EdgeSet::iterator it;
  for (it = p.begin(); it != p.end(); ++it) {
    if ((g[(*it).to_edge_desc(g).first].label < 0) &&
        (g[(*it).in].id == g[v].id)) {
      out += *it;
    }
  }

  return out;
}

EdgeSet flatten(EdgeSets &P, const MyGraph &g) {
  // Concatenate all EdgeSets to a single EdgeSet

  EdgeSet out = EdgeSet(g);
  for (unsigned int i = 0; i < P.size(); ++i) {
    for (unsigned int j = 0; j < P[i].size(); ++j)
      out += P[i][j];
  }
  return out;
}
} // namespace utils
