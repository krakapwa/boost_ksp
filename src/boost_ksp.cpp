#include "boost_ksp.h"

void ksp::add_edge(int n0, int n1, float w=0)
{
    boost::add_edge(n0, n1, Weight(float(w), Label(1)), *G);
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

bool ksp::bellman_ford_shortest_paths(){
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
        weight_map(weight_pmap).
        distance_map(&distance[0]).
        predecessor_map(&parent[0])
        );

    if (r){
        print_distances(distance, *G);
    }
    else
        std::cout << "negative cycle" << std::endl;

    Path shortest = pred_to_path(parent, *G, source_id, sink_id);
    print_path(shortest, *G);
    return r;
}
