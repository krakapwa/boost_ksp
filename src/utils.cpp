#include "utils.h"

VertexPath pred_to_path(std::vector<std::size_t> preds,
                        const MyGraph & g,
                        Vertex source,
                        Vertex sink){
    auto v_index = get(boost::vertex_index, g);

    VertexPath path;
    Vertex current = sink;

    while(g[current].id != g[source].id) {
        path.push_back(current);
        //std::cout << g[current].id << std::endl;
        current=preds[v_index[current]];
    }
    path.push_back(source);

    return path;
}

void print_path(VertexPath path, const MyGraph & g){

    //This prints the path reversed use reverse_iterator and rbegin/rend
    VertexPath::reverse_iterator it;
    //IndexMap index = get(vertex_index, g);
    for (it=path.rbegin(); it != path.rend(); ++it) {

        std::cout << *it << " ";
    }
    std::cout << std::endl;
}

void print_path(EdgePath path, const MyGraph & g){

    for(unsigned int i = 0; i < path.size(); ++i){
        std::cout << "edge: (" << g[source(path[i], g)].name << "," <<
            g[target(path[i], g)].name << ")" << std::endl;
    }

    std::cout << "--------------" << std::endl;
}

void print_all(const MyGraph & g){

    std::cout << "------ edges " << " (n_edges: "
              << num_edges(g) << ") --------" << std::endl;
    EdgeIter ei, ei_end;
    MyGraph::vertex_descriptor v_in;
    MyGraph::vertex_descriptor v_out;

    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
        v_in = source(*ei, g);
        v_out = target(*ei, g);
        std::cout << "edge: (" << g[v_in].id << "," <<
            g[v_out].id << ") / (" << g[v_in].name << "," << g[v_out].name <<
            ") , id: " << g[*ei].id <<
            ", weight: " << g[*ei].weight << std::endl;
    }

    std::cout << "------ vertices " << " (n_vertices: "
              << num_vertices(g) << ") --------" << std::endl;
    std::pair<VertexIter, VertexIter> vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first)
        std::cout << "vertex: " << g[*vp.first].id << ", name: " <<
            g[*vp.first].name << std::endl;
    std::cout << "--------------" << std::endl;

}

EdgePath vertpath_to_edgepath(VertexPath path, const MyGraph & g){
    EdgePath ep;
    std::pair<Edge, bool> e_tmp;

    VertexPath::reverse_iterator it;
    for (it=path.rbegin(); it != path.rend()-1; ++it) {
        //std::cout << "u, v: " << *it << "," << *std::next(it) << std::endl;
        e_tmp = edge(*it, *std::next(it), g);
        ep.push_back(e_tmp.first);

    }
    return ep;
}

void print_dist_pred(std::vector<double> distances,
                     std::vector<Vertex> predecessors,
                     const MyGraph & g){
    std::cout << "distances and parents:\n";
    for (auto v : boost::make_iterator_range(boost::vertices(g))) {
        std::string name = g[v].name;
        std::cout << "distance(" << name << ") = " << distances[v] << ", ";
        std::cout << "parent(" << name << ") = " <<
            g[predecessors[v]].name << "\n";
    }

}
