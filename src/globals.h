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
typedef graph_traits<MyGraph>::edge_descriptor EdgeDesc;
typedef std::vector<graph_traits<MyGraph>::vertex_descriptor> VertexPath;

typedef graph_traits<MyGraph>::vertex_descriptor VertexDesc;
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

// We hold here two vertex descriptors
struct Edge{
  VertexDesc in;
  VertexDesc out;

  Edge(){}

  Edge(const VertexDesc & u,
       const VertexDesc & v){
    in = u;
    out = v;
  }
  
  Edge(const EdgeDesc & e_desc, const MyGraph &a_g){
    in = source(e_desc, a_g);
    out = target(e_desc, a_g);
  }

  std::pair<EdgeDesc, bool> to_edge_desc(const MyGraph & g) const{
    // convert Edge object to EdgeDesc.
    // if graph g doesn't contain such edge, an empty
    // EdgeDesc is returned
    std::pair<EdgeDesc, bool> e = edge(in, out, g);
    if(!e.second){

      BOOST_LOG_TRIVIAL(trace) << "to_edge_desc error."
                              << " edge with vertex descriptors"
                              << " (" << in << "," << out << ")"
                              << "does not exist in graph";
    }
    return e;
  }
};

typedef std::vector<Edge> EdgeVec;

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
    // BOOST_LOG_TRIVIAL(info) << "called copy constructor EdgeSet";
    edges.clear();
    this->g = old.g;
    *this += old;
  }

  EdgeSet &operator+=(const EdgeSet &p) {

    try{
      if(g != p.g) throw 1;
    }
    catch(int e){

      BOOST_LOG_TRIVIAL(trace) << "operator += from EdgeSet to EdgeSet"
                              << " graphs don't match!";

    }
    for(unsigned int i=0; i < p.edges.size(); ++ i){
      *this += p.edges[i];

    }

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

    VertexDesc curr_u;
    VertexDesc curr_v;
    VertexDesc u = e.in;
    VertexDesc v = e.out;
    EdgeVec new_edges;

    for (unsigned int i = 0; i < edges.size(); ++i) {
      // curr_u = source(edges[i], *g);
      curr_u = edges[i].in;
      // curr_v = target(edges[i], *g);
      curr_v = edges[i].out;
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

    return *this;
  }

  void print() const{

    for (unsigned int i = 0; i < edges.size(); ++i) {

      BOOST_LOG_TRIVIAL(trace) << "(" << (*g)[edges[i].in].name << ","
                               << (*g)[edges[i].out].name << ")"
                                << "/"
                               << "(" << (*g)[edges[i].in].id << ","
                               << (*g)[edges[i].out].id << ") ";
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
    VertexDesc curr_u;
    VertexDesc curr_v;

    VertexDesc u = e.in;
    VertexDesc v = e.out;

    for (unsigned int i = 0; i < this->size(); ++i) {
      curr_u = this->edges[i].in;
      curr_v = this->edges[i].out;
      if (((*g)[curr_u].id == (*g)[u].id) && ((*g)[curr_v].id == (*g)[v].id)) {
        this->edges.erase(this->edges.begin() + i);
        break;
      }
    }

    return *this;
  }

  EdgeSet &remove_edge(const VertexDesc &u, const VertexDesc &v) {
    VertexDesc curr_u;
    VertexDesc curr_v;

    for (unsigned int i = 0; i < this->size(); ++i) {
      curr_u = this->edges[i].in;
      curr_v = this->edges[i].out;
      if (((*g)[curr_u].id == (*g)[u].id) && ((*g)[curr_v].id == (*g)[v].id)) {
        this->edges.erase(this->edges.begin() + i);
        break;
      }
    }

    return *this;
  }

  Edge next(const Edge &e) {

    unsigned int ind_e;
    VertexDesc u, v;
    u = e.in;
    v = e.out;
    VertexDesc curr_u, curr_v;

    for (unsigned int i = 0; i < size(); ++i) {
      curr_u = edges[i].in;
      curr_v = edges[i].out;

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

    VertexDesc u, v;
    v = e0.out;
    u = e1.in;
    if ((*g)[u].id == (*g)[v].id)
      return true;
    else
      return false;
  }

  int is_discontinuous() {

    VertexDesc u, v;

    for (unsigned int i = 0; i < size() - 1; ++i) {
      u = edges[i].out;
      v = edges[i + 1].in;
      if ((*g)[u].id != (*g)[v].id)
        return i + 1;
    }
    return -1;
  }

  bool is_contiguous_after(const Edge &e) {

    Edge e_next = next(e);
    VertexDesc u, v;
    v = e.out;
    u = e_next.in;
    if ((*g)[u].id == (*g)[v].id)
      return true;
    else
      return false;
  }

  bool has_vertex(const VertexDesc &u) {

    EdgeSet::iterator it;
    int curr_u_id;
    int curr_v_id;

    for (it = begin(); it != end(); ++it) {
      curr_u_id = ((*g)[(*it).in]).id;
      curr_v_id = ((*g)[(*it).out]).id;
      if ((curr_u_id == (*(this->g))[u].id) ||
          (curr_v_id == (*(this->g))[u].id))
        return true;
    }

    return false;
  }

  bool has_out_vertex(const VertexDesc &u) {

    iterator it;
    vertex_id_type curr_u_id;

    for (it = begin(); it != end(); ++it) {
      curr_u_id = ((*g)[(*it).in]).id;
      if ((curr_u_id == (*(g))[u].id))
        return true;
    }

    return false;
  }

  bool has_edge(const Edge &e) {

    EdgeSet::iterator it;
    vertex_id_type u_id = (*g)[e.in].id;
    vertex_id_type v_id = (*g)[e.out].id;
    vertex_id_type curr_u_id;
    vertex_id_type curr_v_id;

    for (it = begin(); it != end(); ++it) {
      curr_u_id = ((*g)[(*it).in]).id;
      curr_v_id = ((*g)[(*it).out]).id;
      if ((curr_u_id == u_id) && (curr_v_id == v_id))
        return true;
    }
    return false;
  }

  bool sort_descend_labels() {
    // sort edges in descending order according to their label

    // make label vector
    std::vector<int> labels;
    for (unsigned int i = 0; i < edges.size(); ++i)
      labels.push_back((*g)[edges[i].to_edge_desc(*g).first].label);

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
      if ((*g)[edges[i-1].to_edge_desc(*g).first].label >
          (*g)[edges[i].to_edge_desc(*g).first].label)
        return false;

    return true;
  }

  bool has_label(const int &label) {
    // Is there an edge with that label?

    for (unsigned int i = 0; i < edges.size(); ++i){
      if ((*g)[edges[i].to_edge_desc(*g).first].label == label)
        return true;
    }

    return false;
  }

  EdgeSet get_invalid_edges(const MyGraph &g_test) {

    EdgeSet p_invalid(*g);

    for (unsigned int i = 0; i < edges.size(); ++i) {
      std::pair<EdgeDesc, bool> e = edges[i].to_edge_desc(g_test);
      if (!e.second)
        p_invalid += edges[i];
    }

    return p_invalid;
  }

  EdgeSet convert_to_graph(const MyGraph &new_g) {
    // Convert edge descriptors to new graph
    // this works in-place

    EdgeSet p_valid(new_g);

    for (unsigned int i = 0; i < edges.size(); ++i) {
      std::pair<EdgeDesc, bool> e = edges[i].to_edge_desc(new_g);
      if (e.second)
        p_valid += edges[i];
    }

    this->g = p_valid.g;
    this->edges = p_valid.edges;

    return *this;
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

struct EdgeSets {
  std::vector<EdgeSet> sets;
  const MyGraph *g;
  using iterator = std::vector<EdgeSet>::iterator;
  using reverse_iterator = std::vector<EdgeSet>::reverse_iterator;

  EdgeSets(const MyGraph & a_g) { g = &a_g;}

  ~EdgeSets() { sets.clear(); }

  EdgeSets(const EdgeSets &old) {
    // copy constructor
    g = old.g;
    sets = old.sets;
  }

  EdgeSets &operator=(const EdgeSets &other)
  {

    EdgeSets tmp(other);
    this->swap(tmp);
    
    return *this;
  }

  void swap(EdgeSets & other) {
    using std::swap;
    swap(sets, other.sets);
    swap(g, other.g);
  }

  EdgeSets &operator+=(const EdgeSet &new_set) {

    try{
      if(g != new_set.g) throw 1;
    }
    catch(int e){

      BOOST_LOG_TRIVIAL(info) << "operator += from EdgeSet to EdgeSets"
                              << " graphs don't match!";
    }
    sets.push_back(new_set);

    return *this;
  }

  EdgeSets operator+=(const EdgeSets &new_sets) {
    for (unsigned int i = 0; i < new_sets.size(); ++i) 
      *this += new_sets[i];

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

    for (unsigned int i = 0; i < this->sets.size(); ++i) {
      this->sets[i].convert_to_graph(new_g);
    }

    return *this;
  }
};

#endif
