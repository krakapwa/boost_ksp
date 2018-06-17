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
            std::cout << "setting label " << label <<
                " to edge";
            print_edge(es[i], g);
            std::cout << std::endl;
            g[es[i]].label = label;

            std::cout << " edge";
            print_edge(es[i], g);
            std::cout << " label: " << g[es[i]].label << std::endl;
        }
        utils::print_all(g);
    }

    void print_all(const MyGraph & g){

        std::cout << "------ graph: " << g[graph_bundle].name <<
            " --------" << std::endl;
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

    bool edge_is_in_set(Edge e,
                        EdgeSet p,
                        const MyGraph & g_e,
                        const MyGraph & g_p,
                        bool inv_direction=false){

        Vertex v_in_p;
        Vertex v_out_p;

        Vertex v_in_e;
        Vertex v_out_e;

        if(inv_direction){
            v_in_e = target(e, g_e);
            v_out_e = source(e, g_e);
        }
        else{
            v_in_e = source(e, g_e);
            v_out_e = target(e, g_e);
        }

        for(unsigned int i = 0; i < p.size(); ++i){
            v_in_p = source(p[i], g_p);
            v_out_p = target(p[i], g_p);
            if((g_p[v_in_p].id == g_e[v_in_e].id) && (g_p[v_out_p].id == g_e[v_out_e].id))
                return true;
        }

        return false;
    }

    EdgeSet remove_edge_from_set(Vertex v_in,
                                 Vertex v_out,
                                 EdgeSet p,
                                 const MyGraph & g,
                                 bool inv_edge=false){

        Vertex curr_v_in;
        Vertex curr_v_out;

        Vertex v_tmp;

        if(inv_edge){
            v_tmp = v_out;
            v_out = v_in;
            v_in = v_tmp;
        }

        for(unsigned int i = 0; i < p.size(); ++i){
            curr_v_in = source(p[i], g);
            curr_v_out = target(p[i], g);
            if((g[curr_v_in].id == g[v_in].id) &&
               (g[curr_v_out].id == g[v_out].id)){
                //std::cout << "removing edge #" << i << std::endl;
                p.erase(p.begin() + i);
                break;
            }
        }

        return p;
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
                                Vertex v,
                                const MyGraph & g){

        for(unsigned int i=0; i < p.size(); ++i){
            if(g[target(p[i], g)].id == g[v].id)
                return i;
        }

        return -1;
    }

    int find_ind_edge_starting_with(EdgeSet p,
                                    Vertex v,
                                    const MyGraph & g){

        for(unsigned int i=0; i < p.size(); ++i){
            if(g[source(p[i], g)].id == g[v].id)
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
                 Vertex start,
                 Vertex sink,
                 const MyGraph & g,
                 const MyGraph & g_c){

        // Find label of path p (edge is inverted)
        int label_p = g[p[0]].label;

        // Find index of start in p_inter
        EdgeSet to_append;
        int ind_start = find_ind_edge_starting_with(p_inter, start, g);

        Edge curr_edge = p_inter[ind_start];

        // Append edges of p_inter until:
        // (1) sink is reached
        // (2) target of edge is not part of p (based on label)
        // (3) p_inter is "broken"
        print_path(p, g);
        //print_all(g);
        std::cout << "label_stop: " << label_p << std::endl;
        while(true){
            to_append.push_back(edge(source(curr_edge, g_c),
                                     target(curr_edge, g_c), g).first);
            if(g[target(curr_edge, g)].id == g[sink].id)
                break;
            if(g[curr_edge].label == label_p)
                break;
            if(g[p_inter[ind_start]].id_vertex_out !=
               g[p_inter[ind_start+1]].id_vertex_in)
                break;
            ind_start += 1;
            curr_edge = p_inter[ind_start];
        }

        std::cout << "to_append: " << std::endl;
        print_path(to_append, g);
        // Remove edges of to_append from p_inter
        for(unsigned int i=0; i<to_append.size(); ++i)
            p_inter = remove_edge_from_set(to_append[i], p_inter, g);

        // Fetch edges that are "left overs", i.e. interlacing path passed over
        // approach: from vertex target(curr_edge),
        // pass through edges labeled label_p
        // until we find no out_edges with label label_p
        Vertex start_vertex = target(curr_edge, g);
        bool found_a_leftover = true;

        std::cout << "before leftover appending" << std::endl;
        std::cout << "p" << std::endl;
        print_path(p, g);
        print_all(g);
        print_all(g_c);
        std::cout << "---------p_cut-------------" << std::endl;
        print_path(p_cut, g_c);
        print_all(g);

        MyGraph::in_edge_iterator ei, ei_end;
        while(found_a_leftover){
            found_a_leftover = false;
            for (boost::tie(ei, ei_end) = in_edges(start_vertex, g);
                ei != ei_end; ++ei) {

                if((g[*ei].label == label_p) &&
                    !edge_is_in_set(*ei, p_cut, g, g_c, true)){

                    found_a_leftover = true;
                    start_vertex = source(*ei, g);
                    leftovers.push_back(*ei);
                    std::cout << "Leftovers. Appended edge ("
                            << g[source(*ei, g)].name <<  ","
                            << g[target(*ei, g)].name <<  ")"
                            << std::endl;
                }
            }
        }
        //std::cout << "leftovers: " << std::endl;
        //print_path(leftovers, g);

        // append to_append to p

        //std::cout << "to_append: " << std::endl;
        //print_path(to_append, g);
        Vertex append_after_vertex = source(to_append[0], g);
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
            p = remove_edge_from_set(leftovers[i], p, g, false);


        return std::make_tuple(p, p_inter, leftovers, to_append.back());
    }

    EdgeSet translate_edge_set(EdgeSet p,
                               const MyGraph & g_p,
                               const MyGraph & g,
                               bool inv_mode=false){

        EdgeSet p_out;
        Edge e;
        for(unsigned int i=0; i<p.size(); ++i){
            if(inv_mode)
                e = edge(target(p[i], g_p), source(p[i], g_p), g).first;
            else
                e = edge(source(p[i], g_p), target(p[i], g_p), g).first;

            p_out.push_back(e);
        }
        return p_out;
    }

    void invert_edge(Edge e,
                     bool inv_label,
                     bool inv_algebraic_sign,
                     const MyGraph & g_in,
                     MyGraph & g_out){

        Vertex u = source(e, g_in);
        Vertex v = target(e, g_in);
        if(edge(u, v, g_out).second){
            double weight = g_in[e].weight;
            int label = g_in[e].label;
            int id = g_in[e].id;

            if(inv_algebraic_sign)
                weight = -weight;

            if(inv_label)
                label = -label;

            add_edge(g_out,
                    g_in[v].id,
                    g_in[u].id,
                    weight,
                    id,
                    g_in[v].name,
                    g_in[u].name,
                    label);

            boost::remove_edge(u, v, g_out);
        }
    };

    void invert_edges(EdgeSets edge_sets,
                      bool inv_label,
                      bool inv_algebraic_sign,
                      const MyGraph & g_in,
                      MyGraph & g_out){

        EdgeSets::iterator it;
        for (it=edge_sets.begin(); it != edge_sets.end(); ++it) {
            invert_edges(*it, inv_label, inv_algebraic_sign, g_in, g_out);
        }
    }

    void invert_edges(EdgeSet edge_path,
                      bool inv_label,
                      bool inv_algebraic_sign,
                      const MyGraph & g_in,
                      MyGraph & g_out){

        EdgeSet::iterator it;
        for (it=edge_path.begin(); it != edge_path.end(); ++it) {
            invert_edge(*it, inv_label, inv_algebraic_sign, g_in, g_out);
        }
    }

    bool vertex_is_in_set(Vertex v,
                          EdgeSet p,
                          const MyGraph & g){

        EdgeSet::iterator it;
        //int id_v = g[v].id;

        for (it=p.begin(); it != p.end(); ++it) {
            if((g[target(*it, g)].id == g[v].id) || (g[source(*it, g)].id == g[v].id))
                return true;
        }

        return false;

    }

    void set_label_on_path(EdgeSet p,
                            int label,
                            MyGraph & g){

        EdgeSet::iterator it;
        for (it=p.begin(); it != p.end(); ++it) {
            g[*it].label = label;
        }
    }
}
