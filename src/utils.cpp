#include "utils.h"

void print_distances(std::vector<int> distances, const MyGraph & g){
        std::pair<VertexIter, VertexIter> vp;
        for (vp = vertices(g); vp.first != vp.second; ++vp.first)
            std::cout << *vp.first << ": " << distances[*vp.first] << std::endl;
}

VertexPath pred_to_path(std::vector<std::size_t> preds, const MyGraph & g, Vertex source, Vertex sink){
    //std::vector< graph_traits< graph_t >::vertex_descriptor > path;
    VertexPath path;
    Vertex current=sink;

    while(current!=source) {
        path.push_back(current);
        current=preds[current];
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
        std::cout << "edge: " << path[i] << ", id: " << g[path[i]].id << std::endl;
    }

    std::cout << "--------------" << std::endl;
}

void print_all(const MyGraph & g){

    std::cout << "------ edges --------" << std::endl;
    EdgeIter ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        std::cout << "edge id: " << g[*ei].id << ", weight: " << g[*ei].weight << std::endl;

    std::cout << "--------------" << std::endl;

}

EdgePath vertpath_to_edgepath(VertexPath path, const MyGraph & g){
    EdgePath ep;
    std::pair<Edge, bool> e_tmp;

    VertexPath::reverse_iterator it;
    for (it=path.rbegin(); it != path.rend()-1; ++it) {
        std::cout << "u, v: " << *it << "," << *std::next(it) << std::endl;
        e_tmp = edge(*it, *std::next(it), g);
        ep.push_back(e_tmp.first);

    }
    return ep;
}
