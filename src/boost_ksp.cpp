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
    std::cout << "Bellman-Ford" << std::endl;
    std::pair<EdgePath, bool> res;
    res = bellman_ford_shortest_paths();
    print_path(res.first, *G);

    std::cout << "-------------------" << std::endl;
    print_all(*G);
    std::cout << "-------------------" << std::endl;
    std::cout << "Dijkstra" << std::endl;
    res = dijkstra_shortest_paths();
    print_path(res.first, *G);
    return true;
}

std::pair<EdgePath, bool> ksp::dijkstra_shortest_paths(){
    const int nb_vertices = num_vertices(*G);
    std::pair<EdgePath, bool> out;

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
        sink_id,
        weight_map(get(&MyEdge::weight, *G)).
        distance_map(&distance[0]).
        predecessor_map(&parent[0])
        );


    if (r){
        //print_distances(distance, *G);
        VertexPath shortest = pred_to_path(parent, *G, source_id, sink_id);
        //print_path(shortest, *G);
        out.first = vertpath_to_edgepath(shortest, *G);

    }
    else{
        std::cout << "negative cycle" << std::endl;
        out.first = EdgePath();
    }


    out.second = r;
    return out;
}

std::pair<EdgePath, bool> ksp::bellman_ford_shortest_paths(){
    const int nb_vertices = num_vertices(*G);
    std::pair<EdgePath, bool> out;

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
        out.first = vertpath_to_edgepath(shortest, *G);

    }
    else{
        std::cout << "negative cycle" << std::endl;
        out.first = EdgePath();
    }


    out.second = r;
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

void ksp::invert_edges_on_path(EdgePath edge_path, bool inv_algebraic_sign){

    EdgePath::iterator it;
    for (it=edge_path.begin(); it != edge_path.end(); ++it) {
        ksp::invert_edge(*it, inv_algebraic_sign);

    }
}
