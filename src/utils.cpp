#include "utils.h"

void print_distances(std::vector<int> distances, const MyGraph & g){
        std::pair<vertex_iter, vertex_iter> vp;
        for (vp = vertices(g); vp.first != vp.second; ++vp.first)
            std::cout << *vp.first << ": " << distances[*vp.first] << std::endl;
}

Path pred_to_path(std::vector<std::size_t> preds, const MyGraph & g, Vertex source, Vertex sink){
    //std::vector< graph_traits< graph_t >::vertex_descriptor > path;
    Path path;
    Vertex current=sink;

    while(current!=source) {
        path.push_back(current);
        current=preds[current];
    }
    path.push_back(source);

    return path;
}

void print_path(Path path, const MyGraph & g){

    //This prints the path reversed use reverse_iterator and rbegin/rend
    Path::reverse_iterator it;
    IndexMap index = get(vertex_index, g);
    for (it=path.rbegin(); it != path.rend(); ++it) {

        std::cout << index[*it] << " ";
    }
    std::cout << std::endl;
}
