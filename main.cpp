#include <boost/config.hpp>
#include <iostream>
#include <fstream>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

using namespace boost;

enum edge_label_t { edge_label };
enum edge_occupied_t { edge_occupied };
enum vertex_type_t { vertex_type };
enum vertex_frame_t { vertex_frame };
enum vertex_part_t { vertex_part };

namespace boost {
  BOOST_INSTALL_PROPERTY(edge, label);
  BOOST_INSTALL_PROPERTY(edge, occupied);
  BOOST_INSTALL_PROPERTY(vertex, type);
  BOOST_INSTALL_PROPERTY(vertex, frame);
  BOOST_INSTALL_PROPERTY(vertex, part);
}

typedef adjacency_list<
    listS,           // Store out-edges of each vertex in a std::list
    vecS,           // Store vertex set in a std::list
    directedS,    // The file dependency graph is directed
    // vertex properties
        property<vertex_type_t, int,
             property<vertex_frame_t, int,
                        property<vertex_part_t, int> > > ,
    // an edge property
  property<edge_weight_t, float,
           property<edge_label_t, int,
                    property<edge_occupied_t, bool>>>
    > graph_t;
typedef std::pair<int, int> Edge;
typedef graph_traits < graph_t >::vertex_descriptor vertex_t;
typedef graph_traits < graph_t >::edge_descriptor edge_t;
typedef std::vector<vertex_t> vertex_vec;
typedef std::vector<edge_t> path_t;
typedef std::vector<path_t> paths_t;
typedef property_map<graph_t, edge_label_t>::type labels_map_t;
typedef property_map<graph_t, edge_weight_t>::type weights_map_t;
typedef property_map<graph_t, edge_occupied_t>::type occupied_map_t;
typedef property_map<graph_t, vertex_index_t>::type index_map_t;
typename graph_traits < graph_t >::out_edge_iterator ei, ei_end;

void print_path(path_t path, graph_t g){

  weights_map_t weights_map;
  labels_map_t labels_map;
  index_map_t index_map;
  vertex_t v1, v2;
  for(int i = 0; i<path.size();++i){
    v1 = source(path[i],g)

      ;
    v2 = target(path[i],g);
    std::cout << "(" << index_map[v1] << "," << index_map[v2] << ")" <<
      ", w=" << get(weights_map,edge(v1,v2,g).first) << ", l=" << get(labels_map,edge(v1,v2,g).first) << "\t";

  }
  std::cout << std::endl;

}

graph_t interlace(std::vector<path_t> &set, path_t inter, vertex_t sink, graph_t g){

  // Set labels on paths
  bool found;
  edge_t e;
  labels_map_t labels_map = get(edge_label,g);
  occupied_map_t occupied_map = get(edge_occupied,g);
  std::vector<path_t> p_out;
  vertex_t start, end;

  std::vector<path_t> new_set;

  //Label paths
  for(int i = 0; i < set.size(); ++i) {
    path_t this_p = set[i];
    for(int j = 0; j<this_p.size(); ++j){
      start = source(this_p[j],g);
      end = target(this_p[j],g);
      tie(e,found) = edge(start,end,g);
      put(labels_map,e,i+1);
    }
    //print_path(this_p, g);
  }
  std::cout << "done labeling set" << std::endl;

  //Set occupied on edges that are "inverted" (blocked)
  for(int i = 0; i<inter.size(); ++i){
    start = source(inter[i],g);
    end = target(inter[i],g);
    tie(e,found) = edge(start,end,g);
    if(!found){
      tie(e,found) = edge(end,start,g); //invert edge and set occupied flag
      put(occupied_map,e,true);
      //std::cout << "setting occupied edge" << std::endl;
    }
    else {
      std::cout << "labeling interalcing path from " << source(e,g) <<  " to " << target(e,g) << " with label: " << -1 << std::endl;
      put(labels_map,e,-1);
    }
  }

  std::cout << "done setting occupied edges" << std::endl;
  //Interlace paths
  for(int i = 0; i < set.size(); ++i) {
    std::cout << "Interlacing path " << i+1 << std::endl;
    path_t this_p = set[i];
    path_t this_new_p;
    this_new_p.push_back(this_p[0]);

    start = source(this_new_p[0],g);
    end = target(this_new_p[0],g);
    tie(e,found) = edge(start,end,g);
    std::cout << "Adding edge " << source(e,g) <<  " to " << target(e,g) << std::endl;
    this_new_p.push_back(e);
    int label_to_track = labels_map[e];
    bool added_edge = false;

    while(target(e,g) != sink){
      for(tie(ei, ei_end) = out_edges(end, g); ei != ei_end; ++ei) {
        if(labels_map[*ei] == label_to_track){
          if(!occupied_map[*ei]){
            std::cout << "Adding edge " << source(*ei,g) <<  " to " << target(*ei,g) << std::endl;
            this_new_p.push_back(*ei);
            occupied_map[*ei] = true;
            added_edge = true;
          }
        }
      }

      if(!added_edge){
        for(tie(ei, ei_end) = out_edges(end, g); ei != ei_end; ++ei) {
          std::cout << "There is an edge from " << source(*ei,g) <<  " to " << target(*ei,g) << " with label: " << labels_map[*ei] << std::endl;
          if(labels_map[*ei] == -1){
            this_new_p.push_back(*ei);
            std::cout << "Adding edge from " << source(*ei,g) <<  " to " << target(*ei,g) << " with label: " << labels_map[*ei] << std::endl;
            occupied_map[*ei] = true;
            labels_map[*ei] = label_to_track;
          }
        }
      }

      start = source(this_new_p.back(),g);
      end = target(this_new_p.back(),g);
      tie(e,found) = edge(start,end,g);
      added_edge = false;
    }
    new_set.push_back(this_new_p);
  }
  std::cout << "done" << std::endl;

  std::cout << "Interlacing path (pstar) "  << std::endl;
  int label_to_track = -1;
  int new_label = set.size()+1;
  bool added_edge = false;
  path_t this_p = inter;
  path_t new_inter;
  new_inter.push_back(inter[0]);

  start = source(new_inter[0],g);
  end = target(new_inter[0],g);
  tie(e,found) = edge(start,end,g);
  std::cout << "Adding edge " << source(e,g) <<  " to " << target(e,g) << std::endl;
  labels_map[e] = new_label;


  while(target(e,g) != sink){
    for(tie(ei, ei_end) = out_edges(end, g); ei != ei_end; ++ei) {
      if(labels_map[*ei] == label_to_track){
        if(!occupied_map[*ei]){
          std::cout << "Adding edge " << source(*ei,g) <<  " to " << target(*ei,g) << " with label: " << new_label << std::endl;
          new_inter.push_back(*ei);
          occupied_map[*ei] = true;
          labels_map[*ei] = new_label;
          added_edge = true;
        }
      }
    }
    if(!added_edge){
      for(tie(ei, ei_end) = out_edges(end, g); ei != ei_end; ++ei) {
        std::cout << "There is an edge from " << source(*ei,g) <<  " to " << target(*ei,g) << " with label: " << labels_map[*ei] << std::endl;
        if(!occupied_map[*ei]){
          new_inter.push_back(*ei);
          std::cout << "Adding edge from " << source(*ei,g) <<  " to " << target(*ei,g) << " with label: " << new_label << std::endl;
          occupied_map[*ei] = true;
          labels_map[*ei] = new_label;
        }
      }
    }

    start = source(new_inter.back(),g);
    end = target(new_inter.back(),g);
    tie(e,found) = edge(start,end,g);
    added_edge = false;
  }

  new_set.push_back(new_inter);

  set = new_set;

  return g;
}

graph_t cost_transform(graph_t g, std::vector<float> dist){

  weights_map_t weights_map= get(edge_weight,g);

  auto es = edges(g);
  for (auto eit = es.first; eit != es.second; ++eit) {
    float c = weights_map[*eit];
    float s0 = dist[source(*eit,g)];
    float s1 = dist[target(*eit,g)];
    float cp = c + s0 - s1;
    put(weights_map,*eit,cp);
    //std::cout << source(*eit, g) << ' ' << target(*eit, g) << std::endl;
    //std::cout << cp << std::endl;
  }
  return g;
}

graph_t reverseEdgeAndWeightOnPath( path_t path,graph_t g,bool undo_mode=false){
  path_t::iterator it;
  weights_map_t weights_map= get(edge_weight,g);
  float this_weight;
  std::pair<edge_t, bool> new_edge;
  edge_t e;
  bool found;
  vertex_t start, end;

  //std::cout << "lol" << std::endl;
  for(int i = 0; i<path.size(); ++i){
    if(!undo_mode){
      start = source(path[i],g);
      end = target(path[i],g);
    }
    else {
      start = target(path[i],g);
      end = source(path[i],g);
    }
    tie(e,found) = edge(start,end,g);
    this_weight = weights_map[e];
    //std::cout << "(" << start << "," << end << ")" <<
    //  ", w=" << this_weight << "\t";

    remove_edge(start,end,g);
    new_edge = add_edge(end,start,g);
    put(weights_map,new_edge.first,-this_weight);
    //std::cout << "(" << source(new_edge.first,g) << "," << target(new_edge.first,g) << ")" <<
    // ", w=" << weights_map[new_edge.first] << "\t";
  }
  std::cout << std::endl;
  return g;
}

path_t pred_to_path(vertex_vec pred_vec, vertex_t source, vertex_t sink, graph_t g){

  path_t path_out;
  edge_t this_edge;
  vertex_t v = sink;

  for(;;){
    if(v == source) break;
    this_edge = edge(pred_vec[v],v,g).first;
    path_out.push_back(this_edge);
    v = pred_vec[v];
  }

  std::reverse(path_out.begin(),path_out.end());

  return path_out;
}

path_t bellman_shortest_path(vertex_t src, vertex_t snk, graph_t g, std::vector<std::size_t> &parent, std::vector<float> &distance){

  graph_traits < graph_t >::edge_iterator ei, ei_end;
  weights_map_t weights_map = get(edge_weight,g);
  for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    weights_map[*ei] = get(weights_map,*ei);

  for (int i = 0; i < num_vertices(g); ++i)
    parent[i] = i;
  distance[src] = 0;

  bool r = bellman_ford_shortest_paths
    (g, num_vertices(g), weight_map(weights_map).distance_map(&distance[0]).
     predecessor_map(&parent[0]));

  path_t sp = pred_to_path(parent,src,snk,g);

  return sp;
}

int main(int, char *[])
{

  const int num_nodes = 8;
  enum nodes { A, B, C, D, E, F, G, Z };
  char name[] = "ABCDEFGZ";
  typedef std::pair<int,int> Pair;
  Edge edge_array[] = { Edge(A, B), Edge(B, C), Edge(C, D), Edge(D, Z),
                        Edge(A, E), Edge(E, F), Edge(F, Z), Edge(A, G), Edge(G, Z),Edge(E, B),Edge(B, F),Edge(F, D),Edge(C, G)
  };
  float weights[] = { 1, 1, 1, 1, 1, 3, 4, 7, 2, 1, 1, 1, 1};
  int num_edges = sizeof(edge_array) / sizeof(Edge);
  int labels[num_edges] = {0};
  graph_t g(edge_array, edge_array + num_edges, weights, num_nodes);


  vertex_vec p(num_vertices(g));
  std::vector<int> d(num_vertices(g));
  vertex_t snk = vertex(Z, g);
  vertex_t src = vertex(A, g);

  path_t sp;
  std::vector<float> distance(num_vertices(g), (std::numeric_limits < int >::max)());
  std::vector<std::size_t> parent(num_vertices(g));
  sp  = bellman_shortest_path(src, snk, g, parent, distance);
  print_path(sp,g);

  std::cout << "reversing shortest path" << std::endl;
  g = reverseEdgeAndWeightOnPath(sp,g,false);

  g = cost_transform(g, distance);
  std::cout << "p*" << std::endl;

  dijkstra_shortest_paths(g, src,
                            predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, g))).
                            distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));
  path_t inter = pred_to_path(p,src,snk,g);
  print_path(inter,g);

  std::cout << "undo reserving" << std::endl;
  g = reverseEdgeAndWeightOnPath(sp,g,true);
  print_path(sp,g);
  paths_t set;
  set.push_back(sp);
  std::cout << "interlacing:" << std::endl;
  g = interlace(set, inter, snk, g);

  std::cout << "-----" << std::endl;
  print_path(set[0],g);
  std::cout << "-----" << std::endl;
  print_path(set[1],g);

  return EXIT_SUCCESS;
}
