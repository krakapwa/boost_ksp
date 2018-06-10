#include "boost_ksp.h"

void ksp::add_edge(int n0, int n1, float w, int id)
{
    MyEdge e;
    e.weight = w;
    e.label = 1;
    e.id = id;
    boost::add_edge(n0, n1, e, *G);

}

void ksp::set_source_id(int id)
{
    source_id = id;
}

void ksp::set_sink_id(int id)
{
    sink_id = id;
}

void ksp::new_graph(int a_n_vertices){

    std::cout << "Creating graph" << std::endl;
    n_vertices = a_n_vertices;
    G = new MyGraph(n_vertices);

}

shared_ptr<ksp> ksp::create() {

    return shared_ptr<ksp>(new ksp);
}

bool ksp::do_ksp(){

    // Store output of shortest-paths
    std::tuple<EdgePath, bool, std::vector<int>> res;
    EdgePath res_path;
    bool res_ok;
    std::vector<int> res_distance;

    std::cout << "Bellman-Ford" << std::endl;
    res = bellman_ford_shortest_paths();
    std::tie(res_path, res_ok, res_distance) = res;

    print_path(res_path, *G);
    invert_edges(res_path, true);
    print_all(*G);
    cost_transform(res_distance);
    print_all(*G);
    res = dijkstra_shortest_paths();
    std::tie(res_path, res_ok, res_distance) = res;
    print_path(res_path, *G);
    print_all(*G);
    return true;
}

std::tuple<EdgePath, bool, std::vector<int>> ksp::dijkstra_shortest_paths(){
    const int nb_vertices = num_vertices(*G);
    EdgePath out_path;
    std::tuple<EdgePath, bool, std::vector<int>> out;

    // init the distance
    std::vector<int> distance( num_vertices(*G));

    std::vector<std::size_t> predecessors(num_vertices(*G));

    // call to the algorithm
    //bool r = boost::dijkstra_shortest_paths(
    bool r = true;
    boost::dijkstra_shortest_paths(
        *G,
        source_id,
        weight_map(get(&MyEdge::weight, *G)).
        distance_map(make_iterator_property_map(distance.begin(),
                                                       get(vertex_index,*G))).
        predecessor_map(boost::make_iterator_property_map(predecessors.begin(),
                                                          get(vertex_index,*G)))
        );


    if (r){
        //print_distances(distance, *G);
        VertexPath shortest = pred_to_path(predecessors, *G, source_id, sink_id);
        //print_path(shortest, *G);
        out_path = vertpath_to_edgepath(shortest, *G);

    }
    else{
        std::cout << "negative cycle" << std::endl;
        out_path = EdgePath();
    }

    out = std::make_tuple(out_path, r, distance);

    return out;
}

std::tuple<EdgePath, bool, std::vector<int>> ksp::bellman_ford_shortest_paths(){
    std::tuple<EdgePath, bool, std::vector<int>> out;
    EdgePath out_path;

    const int nb_vertices = num_vertices(*G);

    // init the distance
    std::vector<int> distance(nb_vertices, (std::numeric_limits<int>::max)());
    distance[source_id] = 0; // the source is at distance 0

    // init the predecessors (identity function)
    std::vector<std::size_t> parent(nb_vertices);
    for (int i = 0; i < nb_vertices; ++i)
        parent[i] = i;

    // call to the algorithm
    bool r = boost::bellman_ford_shortest_paths(
        *G,
        n_vertices,
        weight_map(get(&MyEdge::weight, *G)).
        distance_map(&distance[0]).
        predecessor_map(&parent[0])
        );

    if (r){
        //print_distances(distance, *G);
        VertexPath shortest = pred_to_path(parent, *G, source_id, sink_id);
        //print_path(shortest, *G);
        out_path = vertpath_to_edgepath(shortest, *G);

    }
    else{
        std::cout << "negative cycle" << std::endl;
        out_path = EdgePath();
    }

    out = std::make_tuple(out_path, r, distance);

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

void ksp::cost_transform(const std::vector<int> & distance){

    EdgeIter ei, ei_end;
    float s_i;
    float s_j;
    float w_ij;
    for (boost::tie(ei, ei_end) = edges(*G); ei != ei_end; ++ei){

        s_i = distance[source(*ei, *G)];
        s_j = distance[target(*ei, *G)];
        w_ij = (*G)[*ei].weight;
        (*G)[*ei].weight = w_ij + s_i - s_j;
    }
}
