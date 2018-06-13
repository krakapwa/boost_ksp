#include "utils.h"

namespace utils{
    Vertex add_vertex(MyGraph & g, int id, std::string str){
        // This passes through all vertices (overhead).

        std::pair<VertexIter, VertexIter> vp;
        for (vp = vertices(g); vp.first != vp.second; ++vp.first)
            if(g[*(vp.first)].id == id){
                //std::cout << "vertex was already there, id: " <<
                //    (*G)[*(vp.first)].id << std::endl;

                return *(vp.first);
            }

        //std::cout << "added new vertex" << std::endl;
        Vertex res = boost::add_vertex(g);
        g[res].name = str;
        g[res].id = id;

        //std::cout << "num_vertices: " << num_vertices(*G) << std::endl;
        return res;
    }

    bool add_edge(MyGraph & g,
                     int n0,
                     int n1,
                     double w,
                     int id=-1,
                     std::string str_0="",
                     std::string str_1="",
                     int label=1){

        // Add two vertices
        Vertex v0 = add_vertex(g, n0, str_0);
        Vertex v1 = add_vertex(g, n1, str_1);

        std::pair<MyGraph::edge_descriptor, bool> e = boost::add_edge(v0, v1, g);

        g[e.first].weight = w;
        g[e.first].label = label;
        g[e.first].id = id;
        g[e.first].id_vertex_in = n0;
        g[e.first].id_vertex_out = n1;

        return true;
    }

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

    void print_edge(Edge e, const MyGraph & g){

        std::cout << "(" << g[source(e, g)].name << "," <<
            g[target(e, g)].name << ") ";
    }
    void print_path(EdgeSet path, const MyGraph & g){

        for(unsigned int i = 0; i < path.size(); ++i){
            print_edge(path[i], g);
        }
        std::cout << std::endl;
    }

    void print_paths(EdgeSets paths, const MyGraph & g){

        for(unsigned int i = 0; i < paths.size(); ++i){
            std::cout << "path " << i << std::endl;
            print_path(paths[i], g);
        }
        std::cout << std::endl;
    }

    EdgeSet get_edges_from_label(const MyGraph & g, int label){
        EdgeIter ei, ei_end;
        EdgeSet out;

        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
            if(g[*ei].label == label)
                out.push_back(*ei);
        }
        return out;
    }

    void set_label_to_all(MyGraph & g, int label){

        EdgeIter ei, ei_end;

        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
            g[*ei].label = label;
        }
    }

    void set_label_to_edges(MyGraph & g, EdgeSet es, int label){

        for(unsigned int i=0; i<es.size(); ++i){
            g[es[i]].label = label;
        }
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
                ", label: " << g[*ei].label <<
                ", weight: " << g[*ei].weight << std::endl;
        }

        std::cout << "------ vertices" << " (n_vertices: "
                << num_vertices(g) << ") --------" << std::endl;
        std::pair<VertexIter, VertexIter> vp;
        for (vp = vertices(g); vp.first != vp.second; ++vp.first)
            std::cout << "vertex: " << g[*vp.first].id << ", name: " <<
                g[*vp.first].name << std::endl;
        std::cout << "--------------" << std::endl;

    }

    EdgeSet vertpath_to_edgepath(VertexPath path, const MyGraph & g){
        EdgeSet ep;
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

    bool edge_is_in_set(Edge e,
                        EdgeSet p,
                        const MyGraph & g,
                        bool inv_direction=false){

        Vertex curr_v_in;
        Vertex curr_v_out;

        Vertex v_in;
        Vertex v_out;

        if(inv_direction){
            v_in = target(e, g);
            v_out = source(e, g);
        }
        else{
            v_in = source(e, g);
            v_out = target(e, g);
        }

        for(unsigned int i = 0; i < p.size(); ++i){
            curr_v_in = source(p[i], g);
            curr_v_out = target(p[i], g);
            if((g[curr_v_in].id == g[v_in].id) && (g[curr_v_out].id == g[v_out].id))
                return true;
        }

        return false;
    }

    EdgeSet remove_edge_from_set(Edge e,
                                EdgeSet p,
                                const MyGraph & g,
                                bool inv_edge=false){

        Vertex curr_v_in;
        Vertex curr_v_out;

        Vertex v_in;
        Vertex v_out;

        if(inv_edge){
            v_in = target(e, g);
            v_out = source(e, g);
        }
        else{
            v_in = source(e, g);
            v_out = target(e, g);
        }

        for(unsigned int i = 0; i < p.size(); ++i){
            curr_v_in = source(p[i], g);
            curr_v_out = target(p[i], g);
            if((g[curr_v_in].id == g[v_in].id) && (g[curr_v_out].id == g[v_out].id)){
                //std::cout << "removing edge #" << i << std::endl;
                p.erase(p.begin() + i);
                break;
            }
        }

        return p;
    }
    int find_ind_edge_ending_with(EdgeSet p,
                                MyVertex v,
                                const MyGraph & g){

        for(unsigned int i=0; i < p.size(); ++i){
            if(g[target(p[i], g)].id == v.id)
                return i;
        }

        return -1;
    }

    int find_ind_edge_starting_with(EdgeSet p,
                                    MyVertex v,
                                    const MyGraph & g){

        for(unsigned int i=0; i < p.size(); ++i){
            if(g[source(p[i], g)].id == v.id)
                return i;
        }

        //std::cout << "find_ind_edge_starting_with: Could find requested edge.";
        //std::cout << "Something bad happened" << std::endl;
        return -1;
    }

    std::tuple<EdgeSet, EdgeSet, EdgeSet, Edge>
    append_inter(EdgeSet p,
                EdgeSet p_inter,
                EdgeSet p_cut,
                EdgeSet leftovers,
                MyVertex start,
                MyVertex sink,
                const MyGraph & g){

        // Find index of start in p_inter
        EdgeSet to_append;
        int ind_start = find_ind_edge_starting_with(p_inter, start, g);

        Edge curr_edge = p_inter[ind_start];

        // Append edges of p_inter until: (1) sink is not reached
        // (2) target of edge is not part of p (case where we stop)
        while(true){
            to_append.push_back(curr_edge);
            if(g[target(curr_edge, g)].id == sink.id)
                break;
            if(find_ind_edge_starting_with(p, g[target(curr_edge, g)], g) != -1)
                break;
            ind_start += 1;
            curr_edge = p_inter[ind_start];
        }

        // Remove edges of to_append from p_inter
        for(unsigned int i=0; i<to_append.size(); ++i)
            p_inter = remove_edge_from_set(to_append[i], p_inter, g);

        // Fetch edges that are "left overs", i.e. interlacing path passed over
        // will be added to p_inter
        // approach: from vertex target(curr_edge), pass through edges labeled -1
        // until we find no out_edges with label -1
        Vertex start_vertex = target(curr_edge, g);
        MyGraph::out_edge_iterator ei, ei_end;
        bool found_a_leftover = true;

        while(found_a_leftover){
            found_a_leftover = false;
            for (boost::tie(ei, ei_end) = out_edges(start_vertex, g);
                ei != ei_end; ++ei) {
                if(g[*ei].label == -1){
                    if(edge_is_in_set(*ei, p_cut, g, false)){
                        //std::cout << "hit a cut edge, breaking" << std::endl;
                        break;
                    }
                    found_a_leftover = true;
                    start_vertex = target(*ei, g);
                    leftovers.push_back(*ei);
                    //std::cout << "Leftovers. Appended edge from ("
                    //          << g[source(*ei, g)].name <<  ","
                    //          << g[target(*ei, g)].name <<  ")"
                    //          << std::endl;
                }
            }
        }
        //std::cout << "leftovers: " << std::endl;
        //print_path(leftovers, g);

        // append to_append to p

        //std::cout << "to_append: " << std::endl;
        //print_path(to_append, g);
        MyVertex append_after_vertex = g[source(to_append[0], g)];
        //std::cout << "append_after_vertex: " << append_after_vertex.name << std::endl;
        //std::cout << "p before appending to_append: " << std::endl;
        //print_path(p, g);
        int ind_to_append = find_ind_edge_ending_with(p,
                                                        append_after_vertex,
                                                        g);
        ind_to_append += 1;
        //std::cout << "ind_to_append: " << ind_to_append << std::endl;

        for(unsigned int i=0; i<to_append.size(); ++i)
            p.insert(p.begin() + ind_to_append + i, to_append[i]);
        //std::cout << "p after appending to_append: " << std::endl;
        //print_path(p, g);

        // remove left overs from p
        for(unsigned int i=0; i<leftovers.size(); ++i)
            p = remove_edge_from_set(leftovers[i], p, g, true);


        return std::make_tuple(p, p_inter, leftovers, to_append.back());
    }


    void invert_edge(Edge e,
                        bool inv_label,
                        bool inv_algebraic_sign,
                        MyGraph & g){

        Vertex u = source(e, g);
        Vertex v = target(e, g);
        double weight = g[e].weight;
        int label = g[e].label;
        int id = g[e].id;

        if(inv_algebraic_sign)
            weight = -weight;

        if(inv_label)
            label = -label;

        add_edge(g,
                 g[v].id,
                 g[u].id,
                 weight,
                 id,
                 g[v].name,
                 g[u].name,
                 label);

        boost::remove_edge(u, v, g);
    };

    void invert_edges(EdgeSet edge_path,
                    bool inv_label,
                    bool inv_algebraic_sign,
                    MyGraph & g){

        EdgeSet::iterator it;
        for (it=edge_path.begin(); it != edge_path.end(); ++it) {
            invert_edge(*it, inv_label, inv_algebraic_sign, g);
        }
    }

    void invert_edges(EdgeSets set,
                    bool inv_label,
                    bool inv_algebraic_sign,
                    MyGraph & g){

        for(unsigned int i=0; i<set.size(); ++i){
            invert_edges(set[i], inv_label, inv_algebraic_sign, g);
        }
    }
}
