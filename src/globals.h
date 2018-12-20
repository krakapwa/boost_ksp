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

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/property_map/transform_value_property_map.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/shared_ptr.hpp>
#include <numeric>
#include <tuple>
#include <iostream>

namespace bp = boost::python;
namespace bn = boost::python::numpy;
namespace bl = boost::log;

using namespace boost;

enum logSeverityLevel {
  trace,   // 5
  debug,   // 4
  info,    // 3
  warning, // 2
  error    // 1
};

enum edge_myweight_t { edge_myweight };
enum edge_label_t { edge_label };
enum edge_id_t { edge_id };

typedef long int vertex_id_type;
typedef long int edge_id_type;

struct GraphProperty {
  std::string name;
};

struct EdgeProperty {
  float weight;
  int label;       // used during optimization
  edge_id_type id; // will return solution as list of ids, make them unique!
  vertex_id_type id_vertex_in;
  vertex_id_type id_vertex_out; // ids of input and output vertices
};

struct VertexProperty {
  std::string name;
  vertex_id_type id;

  bool operator==(const VertexProperty & v) {
    if (id == v.id)
      return true;
    else
      return false;
  }

  bool operator!=(const VertexProperty & v) {
    if (id != v.id)
      return true;
    else
      return false;
  }
};

typedef adjacency_list<vecS, vecS, directedS, VertexProperty, EdgeProperty,
                       GraphProperty>
    MyGraph;
// typedef adjacency_list<vecS, vecS, bidirectionalS,
//     VertexProperty, EdgeProperty, GraphProperty> MyGraph;

typedef graph_traits<MyGraph>::vertex_iterator VertexIter;
typedef graph_traits<MyGraph>::edge_descriptor Edge;
typedef std::vector<graph_traits<MyGraph>::vertex_descriptor> VertexPath;

typedef graph_traits<MyGraph>::vertex_descriptor Vertex;
typedef std::vector<graph_traits<MyGraph>::edge_descriptor> EdgeVec;
typedef graph_traits<MyGraph>::edge_iterator EdgeIter;
typedef graph_traits<MyGraph>::out_edge_iterator OutEdgeIter;

template <typename T>
std::vector<std::size_t> idx_desc_sort(const std::vector<T> &v) {
  std::vector<std::size_t> result(v.size());
  std::iota(std::begin(result), std::end(result), 0);
  std::sort(std::begin(result), std::end(result),
            [&v](const auto &lhs, const auto &rhs) { return v[lhs] > v[rhs]; });
  return result;
}

// We'll make extensive use of edge sets
struct EdgeSet {
  EdgeVec edges;
  const MyGraph *g;
  using iterator = EdgeVec::iterator;
  using reverse_iterator = EdgeVec::reverse_iterator;

  // Default constructor. Must assign a graph!
  EdgeSet(const MyGraph &a_g) { g = &a_g; }

  ~EdgeSet() { edges.clear(); }

  // copy constructor
  EdgeSet(const EdgeSet &old) {
    BOOST_LOG_TRIVIAL(info) << "called copy constructor EdgeSet";
    edges.clear();
    this->g = old.g;
    *this += old;
  }

  EdgeSet &operator+=(const EdgeSet &p) {

    BOOST_LOG_TRIVIAL(debug) << "in append EdgeSet -> EdgeSet";
    // BOOST_LOG_TRIVIAL(info) << "this->g: " << g
    //                         << " p.g: " << p.g;
    try{
      if(g != p.g) throw 1;
    }
    catch(int e){

      BOOST_LOG_TRIVIAL(info) << "exception " << e;
      BOOST_LOG_TRIVIAL(info) << "operator += from EdgeSet to EdgeSet"
                              << " graphs don't match!";

    }
    BOOST_LOG_TRIVIAL(info) << "p.edges.size(): " << p.edges.size();
    BOOST_LOG_TRIVIAL(info) << "edges.size(): " << edges.size();
    // edges.insert(this->edges.end(), p.edges.begin(), p.edges.end());
    for(unsigned int i=0; i < p.edges.size(); ++ i){
      // BOOST_LOG_TRIVIAL(info) << "p.edges[" << i << "]: " << p.edges[i];
      *this += p.edges[i];

    }
    BOOST_LOG_TRIVIAL(debug) << "inserted!";

    return *this;
  }

  EdgeSet &operator+=(const Edge &e) {

    edges.push_back(e);

    return *this;
  }

  EdgeSet operator-=(const EdgeSet & p) {

    for (unsigned int i = 0; i < p.size(); ++i) {
      *this -= p[i];
    }

    return *this;
  }

  EdgeSet operator-=(const Edge &e) {

    Vertex curr_u;
    Vertex curr_v;
    Vertex u = source(e, *g);
    Vertex v = target(e, *g);
    EdgeVec new_edges;

    for (unsigned int i = 0; i < edges.size(); ++i) {
      curr_u = source(edges[i], *g);
      curr_v = target(edges[i], *g);
      if (!(((*g)[u].id == (*g)[curr_u].id) &&
            ((*g)[v].id == (*g)[curr_v].id))) {
        new_edges.push_back(edges[i]);
      } else {
        BOOST_LOG_TRIVIAL(trace) << "removed 1 edge";
      }
    }

    edges = new_edges;

    return *this;
  }

  const Edge &operator[](const edge_id_type &i) const { return edges[i]; }

  EdgeSet & operator=(const EdgeSet &new_set) {
    g = new_set.g;
    edges.clear();
    edges = new_set.edges;
    BOOST_LOG_TRIVIAL(trace) << "EdgeSet operator= returning";

    return *this;
  }

  void print() const{

    for (unsigned int i = 0; i < edges.size(); ++i) {

      BOOST_LOG_TRIVIAL(trace) << "(" << (*g)[source(edges[i], *g)].name << ","
                               << (*g)[target(edges[i], *g)].name << ")"
                                << "/"
                               << "(" << (*g)[source(edges[i], *g)].id << ","
                               << (*g)[target(edges[i], *g)].id << ") "
                               << "label: " << (*g)[edges[i]].label
                               << " weight: " << (*g)[edges[i]].weight;
    }
    
  }

  unsigned int size() const { return edges.size(); }

  Edge back() { return edges.back(); }

  iterator begin() { return edges.begin(); }

  iterator end() { return edges.end(); }

  reverse_iterator rbegin() { return edges.rbegin(); }

  reverse_iterator rend() { return edges.rend(); }

  EdgeSet &insert(EdgeSet &p, const int &ind) {

    for (unsigned int i = 0; i < p.size(); ++i)
      edges.insert(edges.begin() + ind + i, p[i]);

    return *this;
  }

  EdgeSet &remove_edge(const Edge &e) {
    Vertex curr_u;
    Vertex curr_v;

    Vertex u = source(e, *g);
    Vertex v = target(e, *g);

    for (unsigned int i = 0; i < this->size(); ++i) {
      curr_u = source(this->edges[i], *g);
      curr_v = target(this->edges[i], *g);
      if (((*g)[curr_u].id == (*g)[u].id) && ((*g)[curr_v].id == (*g)[v].id)) {
        this->edges.erase(this->edges.begin() + i);
        break;
      }
    }

    return *this;
  }

  EdgeSet &remove_edge(const Vertex &u, const Vertex &v) {
    Vertex curr_u;
    Vertex curr_v;

    for (unsigned int i = 0; i < this->size(); ++i) {
      curr_u = source(this->edges[i], *g);
      curr_v = target(this->edges[i], *g);
      if (((*g)[curr_u].id == (*g)[u].id) && ((*g)[curr_v].id == (*g)[v].id)) {
        this->edges.erase(this->edges.begin() + i);
        break;
      }
    }

    return *this;
  }

  Edge next(const Edge &e) {

    unsigned int ind_e;
    Vertex u, v;
    u = source(e, *g);
    v = target(e, *g);
    Vertex curr_u, curr_v;

    for (unsigned int i = 0; i < size(); ++i) {
      curr_u = source(edges[i], *g);
      curr_v = target(edges[i], *g);

      if ((*g)[curr_u].id == (*g)[u].id && (*g)[curr_v].id == (*g)[v].id) {
        ind_e = i;
        break;
      }
    }

    if (ind_e < size() - 1) // we got last edge, return first edge
      return edges[ind_e + 1];
    else
      return edges[0];
  }

  bool are_contiguous(const Edge &e0, const Edge &e1) {

    Vertex u, v;
    v = target(e0, *g);
    u = source(e1, *g);
    if ((*g)[u].id == (*g)[v].id)
      return true;
    else
      return false;
  }

  int is_discontinuous() {

    Vertex u, v;

    for (unsigned int i = 0; i < size() - 1; ++i) {
      u = target(edges[i], *g);
      v = source(edges[i + 1], *g);
      if ((*g)[u].id != (*g)[v].id)
        return i + 1;
    }
    return -1;
  }

  bool is_contiguous_after(const Edge &e) {

    Edge e_next = next(e);
    Vertex u, v;
    v = target(e, *g);
    u = source(e_next, *g);
    if ((*g)[u].id == (*g)[v].id)
      return true;
    else
      return false;
  }

  bool has_vertex(const Vertex &u) {

    EdgeSet::iterator it;
    int curr_u_id;
    int curr_v_id;

    for (it = begin(); it != end(); ++it) {
      curr_u_id = ((*g)[source(*it, *g)]).id;
      curr_v_id = ((*g)[target(*it, *g)]).id;
      if ((curr_u_id == (*(this->g))[u].id) ||
          (curr_v_id == (*(this->g))[u].id))
        return true;
    }

    return false;
  }

  bool has_out_vertex(const Vertex &u) {

    iterator it;
    vertex_id_type curr_u_id;

    for (it = begin(); it != end(); ++it) {
      curr_u_id = ((*g)[source(*it, *g)]).id;
      if ((curr_u_id == (*(g))[u].id))
        return true;
    }

    return false;
  }

  bool has_edge(const Edge &e) {

    EdgeSet::iterator it;
    vertex_id_type u_id = (*g)[source(e, *g)].id;
    vertex_id_type v_id = (*g)[target(e, *g)].id;
    vertex_id_type curr_u_id;
    vertex_id_type curr_v_id;

    for (it = begin(); it != end(); ++it) {
      curr_u_id = ((*g)[source(*it, *g)]).id;
      curr_v_id = ((*g)[target(*it, *g)]).id;
      if ((curr_u_id == u_id) && (curr_v_id == v_id))
        return true;
    }
    return false;
  }

  EdgeSet &remove_label(const int &label) {

    EdgeSet out(*g);
    for (unsigned int i = 0; i < edges.size(); ++i)
      if ((*g)[edges[i]].label != label)
        out += edges[i];

    edges = out.edges;

    return *this;
  }

  EdgeSet &keep_label(const int &label) {

    EdgeSet out(*g);
    for (unsigned int i = 0; i < edges.size(); ++i)
      if ((*g)[edges[i]].label == label)
        out += edges[i];

    edges = out.edges;

    return *this;
  }

  bool sort_descend_labels() {
    // sort edges in descending order according to their label

    // make label vector
    std::vector<int> labels;
    for (unsigned int i = 0; i < edges.size(); ++i)
      labels.push_back((*g)[edges[i]].label);

    // get indices that sort labels in descending order
    auto idxs = idx_desc_sort(labels);

    EdgeVec new_edges;
    for (unsigned int i = 0; i < idxs.size(); ++i)
      new_edges.push_back(edges[idxs[i]]);

    edges = new_edges;

    return true;
  }

  bool are_label_sorted() {
    // are edges sorted in descending order according to their label?

    for (unsigned int i = 1; i < edges.size(); ++i)
      if ((*g)[edges[i]].label > (*g)[edges[i - 1]].label)
        return false;

    return true;
  }

  bool has_label(const int &label) {
    // Is there an edge with that label?

    for (unsigned int i = 0; i < edges.size(); ++i)
      if ((*g)[edges[i]].label == label)
        return true;

    return false;
  }

  EdgeSet get_invalid_edges(const MyGraph &g_test) {

    EdgeSet p_invalid(*g);
    std::pair<Edge, bool> e;

    Vertex u, v;
    for (unsigned int i = 0; i < edges.size(); ++i) {
      u = source(edges[i], *g);
      v = target(edges[i], *g);
      e = edge(u, v, g_test);
      if (!e.second)
        p_invalid += edges[i];
    }

    return p_invalid;
  }

  EdgeSet convert_to_graph(const MyGraph &new_g) {
    // Convert edge descriptors to new graph
    // this works in-place

    EdgeSet p_valid(new_g);
    std::pair<Edge, bool> e;

    Vertex u, v;

    for (unsigned int i = 0; i < edges.size(); ++i) {
      u = source(edges[i], *g);
      v = target(edges[i], *g);
      e = edge(u, v, new_g);
      if (e.second)
        p_valid += e.first;
    }

    this->g = p_valid.g;
    this->edges = p_valid.edges;

    return *this;
  }

  std::pair<Edge, bool> convert_edge(Edge e, const MyGraph &g_p,
                                     bool inv_mode) {
    // Convert edge to new graph
    // In place conversion

    std::pair<Edge, bool> e_out;
    if (inv_mode)
      e_out = edge(target(e, g_p), source(e, g_p), *g);
    else
      e_out = edge(source(e, g_p), target(e, g_p), *g);

    return e_out;
  }
};

inline EdgeSet operator+(const EdgeSet &p0, const EdgeSet &p1) {
  // Concatenate two sets

  EdgeSet p_out(*(p0.g));
  p_out += p0;
  p_out += p1;

  return p_out;
}

inline EdgeSet operator+(const EdgeSet &p0, const Edge &e) {
  // Concatenate edge e to set p0

  EdgeSet p_out(*(p0.g));
  p_out += p0;
  p_out += e;

  return p_out;
}

inline EdgeSet operator-(const EdgeSet & p0, const EdgeSet & p1) {
  // Set substraction

  EdgeSet p_out(*(p0.g));
  p_out += p0;
  for (unsigned int i = 0; i < p1.size(); ++i) {
    if (p_out.has_edge(p1[i]))
      p_out -= p1[i];
  }

  return p_out;
}

class LabelSorter {
  MyGraph g;

public:
  LabelSorter(MyGraph a_g) { g = a_g; }
  bool compareLabels(const Edge e0, const Edge e1, MyGraph g) const {

    BOOST_LOG_TRIVIAL(trace)
        << "compareLabels (e0, e1): " << g[e0].label << "," << g[e1].label
        << " return " << (g[e0].label > g[e1].label);
    return g[e0].label > g[e1].label;
  }
  bool operator()(const Edge & e0, const Edge & e1) const {
    return compareLabels(e0, e1, g);
  }
};

struct EdgeSets {
  std::vector<EdgeSet> sets;
  const MyGraph *g;
  using iterator = std::vector<EdgeSet>::iterator;
  using reverse_iterator = std::vector<EdgeSet>::reverse_iterator;

  EdgeSets(const MyGraph & a_g) { g = &a_g;}

  ~EdgeSets() { sets.clear(); }

  EdgeSets(const EdgeSets &old) {
    // copy constructor
    BOOST_LOG_TRIVIAL(debug) << "in copy constructor EdgeSets";
    g = old.g;
    sets = old.sets;
    // BOOST_LOG_TRIVIAL(trace) << "in copy constructor. done append";
  }

  EdgeSets &operator=(const EdgeSets &other)
  {
    BOOST_LOG_TRIVIAL(debug) << "in operator=";

    EdgeSets tmp(other);
    this->swap(tmp);
    
    BOOST_LOG_TRIVIAL(debug) << "copied";
    return *this;
  }

  void swap(EdgeSets & other) {
    using std::swap;
    swap(sets, other.sets);
    swap(g, other.g);
  }

  EdgeSets &operator+=(const EdgeSet &new_set) {

    BOOST_LOG_TRIVIAL(debug) << "in append EdgeSet -> EdgeSets";
      // BOOST_LOG_TRIVIAL(info) << "this->g: " << g
      //                         << " new_set.g: " << new_set.g;

    try{
      if(g != new_set.g) throw 1;
    }
    catch(int e){

      BOOST_LOG_TRIVIAL(info) << "exception " << e;
      BOOST_LOG_TRIVIAL(info) << "operator += from EdgeSet to EdgeSets"
                              << " graphs don't match!";
    }
    BOOST_LOG_TRIVIAL(debug) << "will push_back EdgeSet to EdgeSets";
    sets.push_back(new_set);
    BOOST_LOG_TRIVIAL(debug) << "done append EdgeSet -> EdgeSets";

    return *this;
  }

  EdgeSets operator+=(const EdgeSets &new_sets) {
    BOOST_LOG_TRIVIAL(debug) << "in append EdgeSets -> EdgeSets";
    BOOST_LOG_TRIVIAL(trace) << "new_sets.size():  " << new_sets.size();
    BOOST_LOG_TRIVIAL(trace) << "this->size(): " << this->size();
    // new_sets.print();

    for (unsigned int i = 0; i < new_sets.size(); ++i) {
      BOOST_LOG_TRIVIAL(trace) << "in append loop";
      *this += new_sets[i];
    }
    return *this;
  }

  void print() const{
    for (unsigned int i = 0; i < sets.size(); ++i) {
      sets[i].print();
    }

  }

  void insert(const EdgeSet &new_set) {
    // push new EdgeSet at beginning
    sets.insert(sets.begin(), new_set);
  }

  unsigned int size() const { return sets.size(); }

  const EdgeSet & operator[](const int & i) const {

    if(i >= sets.size()){
      BOOST_LOG_TRIVIAL(debug) << "operator[]: asking for invalid index: "
                               << i << ". size is: " << sets.size();
}

    return sets[i];
  }

  EdgeSet back() { return sets.back(); }

  iterator begin() { return sets.begin(); }

  iterator end() { return sets.end(); }

  reverse_iterator rbegin() { return sets.rbegin(); }

  reverse_iterator rend() { return sets.rend(); }

  EdgeSets convert_to_graph(const MyGraph &new_g) {
    // This converts in-place

    BOOST_LOG_TRIVIAL(debug) << "EdgeSets convert_to_graph";

    for (unsigned int i = 0; i < this->sets.size(); ++i) {
      this->sets[i].convert_to_graph(new_g);
    }

    return *this;
  }
};

#endif
