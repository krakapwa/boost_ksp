#include "boost_ksp.h"

Vertex ksp::my_add_vertex(int id, std::string str){
    // This passes through all vertices (overhead).

    std::pair<VertexIter, VertexIter> vp;
    for (vp = vertices(*G); vp.first != vp.second; ++vp.first)
        if((*G)[*(vp.first)].id == id){
            //std::cout << "vertex was already there, id: " <<
            //    (*G)[*(vp.first)].id << std::endl;

            return *(vp.first);
        }

    //std::cout << "added new vertex" << std::endl;
    Vertex res = boost::add_vertex(*G);
    (*G)[res].name = str;
    (*G)[res].id = id;

    //std::cout << "num_vertices: " << num_vertices(*G) << std::endl;
    return res;
}

void ksp::add_edge(int n0,
                   int n1,
                   double w,
                   int id,
                   std::string str_0="",
                   std::string str_1="")
{

    // Add two vertices
    Vertex v0 = my_add_vertex(n0, str_0);
    Vertex v1 = my_add_vertex(n1, str_1);

    std::pair<MyGraph::edge_descriptor, bool> e = boost::add_edge(v0, v1, *G);

    (*G)[e.first].weight = w;
    (*G)[e.first].label = 1;
    (*G)[e.first].id = id;

}

void ksp::set_source(int id, std::string str)
{
    source_vertex = my_add_vertex(id, str);

}

void ksp::set_sink(int id, std::string str)
{
    sink_vertex = my_add_vertex(id, str);
}

void ksp::new_graph(int a_n_vertices=0){

    n_vertices = a_n_vertices;
    std::cout << "Creating graph with n_vertices: " <<
              n_vertices << std::endl;
    G = new MyGraph(n_vertices);

}

shared_ptr<ksp> ksp::create() {

    return shared_ptr<ksp>(new ksp);
}

bool ksp::do_ksp(){

    // Store output of shortest-paths
    ShortestPathRes res;
    EdgePath res_path;
    bool res_ok;
    std::vector<double> res_distance;


    std::cout << "Bellman-Ford" << std::endl;
    res = bellman_ford_shortest_paths();
    std::tie(res_path, res_ok, res_distance) = res;
    print_path(res_path, *G);

    invert_edges(res_path, true);
    print_all(*G);
    cost_transform(res_distance);
    std::cout << "Dijkstra" << std::endl;
    res = dijkstra_shortest_paths();
    std::tie(res_path, res_ok, res_distance) = res;
    print_path(res_path, *G);
    return true;
}

ShortestPathRes ksp::dijkstra_shortest_paths(){

    EdgePath out_path;
    ShortestPathRes out;

    // init the distance
    std::vector<double> distances( num_vertices(*G));

    std::vector<std::size_t> predecessors(num_vertices(*G));

    // call to the algorithm
    //bool r = boost::dijkstra_shortest_paths(
    bool r = true;
    boost::dijkstra_shortest_paths(
        *G,
        source_vertex,
        weight_map(get(&MyEdge::weight, *G)).
        distance_map(make_iterator_property_map(distances.begin(),
                                                       get(vertex_index,*G))).
        predecessor_map(boost::make_iterator_property_map(predecessors.begin(),
                                                          get(vertex_index,*G)))
        );


    if (r){
        //print_distances(distance, *G);
        VertexPath shortest = pred_to_path(predecessors, *G,
                                           source_vertex, sink_vertex);
        //print_path(shortest, *G);
        out_path = vertpath_to_edgepath(shortest, *G);

    }
    else{
        std::cout << "negative cycle" << std::endl;
        out_path = EdgePath();
    }

    out = std::make_tuple(out_path, r, distances);

    return out;
}

ShortestPathRes  ksp::bellman_ford_shortest_paths(){
    ShortestPathRes   out;
    EdgePath out_path;
    std::pair<VertexIter, VertexIter> vp;
    auto v_index = get(boost::vertex_index, *G);
    auto weight = boost::make_transform_value_property_map(
        std::mem_fn(&MyEdge::weight),
        get(boost::edge_bundle, *G));

    // init
    std::vector<Vertex> predecessors(num_vertices(*G));
    //for (vp = vertices(*G); vp.first != vp.second; ++vp.first)
    //    predecessors[v_index[*vp.first]] = v_index[*vp.first];

    std::vector<double> distances(num_vertices(*G),
                                  (std::numeric_limits<double>::max)());
    distances[v_index[source_vertex]] = 0;

    bool r = boost::bellman_ford_shortest_paths(
        *G,
        num_vertices(*G),
        weight_map(weight).
        distance_map(make_iterator_property_map(distances.begin(),
                                                v_index)).
        predecessor_map(make_iterator_property_map(predecessors.begin(),
                                                   v_index))
        );

    if (r){
        //print_dist_pred(distances, predecessors, *G);

        VertexPath shortest = pred_to_path(predecessors, *G,
                                           source_vertex, sink_vertex);
        out_path = vertpath_to_edgepath(shortest, *G);

    }
    else{
        std::cout << "negative cycle" << std::endl;
        out_path = EdgePath();
    }

    out = std::make_tuple(out_path, r, distances);

    return out;
}

void ksp::print_edges(){

    print_all(*G);
};

void ksp::invert_edge(Edge e, bool inv_algebraic_sign){

    MyEdge e_new;
    e_new.label = (*G)[e].label;
    e_new.id = (*G)[e].id;
    int u = source(e, *G);
    int v = target(e, *G);

    if(inv_algebraic_sign)
        e_new.weight = -(*G)[e].weight;


    boost::add_edge(v, u, e_new, *G);
    boost::remove_edge(u, v, *G);

};

void ksp::invert_edges(EdgePath edge_path, bool inv_algebraic_sign){

    EdgePath::iterator it;
    for (it=edge_path.begin(); it != edge_path.end(); ++it) {
        ksp::invert_edge(*it, inv_algebraic_sign);

    }
}

void ksp::cost_transform(const std::vector<double> & distance){

    EdgeIter ei, ei_end;
    double s_i;
    double s_j;
    double w_ij;
    for (boost::tie(ei, ei_end) = edges(*G); ei != ei_end; ++ei){

        s_i = distance[source(*ei, *G)];
        s_j = distance[target(*ei, *G)];
        w_ij = (*G)[*ei].weight;
        (*G)[*ei].weight = w_ij + s_i - s_j;
    }
}
