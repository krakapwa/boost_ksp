/*
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributor(s): Laurent Lejeune (laurent.lejeune@artorg.unibe.ch).
 *
 */

#include "utils.h"

namespace utils{
    // This passes through all vertices before adding (overhead).
    Vertex add_vertex(MyGraph & g, int id, std::string str){


        std::pair<VertexIter, VertexIter> vp;
        for (vp = vertices(g); vp.first != vp.second; ++vp.first)
            if(g[*(vp.first)].id == id){
                //vertex was already there
                return *(vp.first);
            }

        // addnew vertex
        Vertex res = boost::add_vertex(g);
        g[res].name = str;
        g[res].id = id;

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

    // Converts from predecessor map to vector of edge descriptors
    std::tuple<VertexPath,bool> pred_to_path(std::vector<std::size_t> preds,
                            const MyGraph & g,
                            Vertex source,
                            Vertex sink){

        auto v_index = get(boost::vertex_index, g);

        VertexPath path;
        Vertex current = sink;
        Vertex last;
        bool ok = true;

        path.push_back(current);

        while(g[current].id != g[source].id) {

            last = current;
            current=preds[v_index[current]];

            if(g[current].id == g[last].id){
                ok = false;
                break;
            }

            path.push_back(current);
        }

        //path.push_back(source);

        return std::make_tuple(path, ok);
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
        BOOST_LOG_TRIVIAL(info) << "(" << g[source(e, g)].name
                                 << ","
                                 << g[target(e, g)].name << ") ";
    }

    void print_path(EdgeSet path, const MyGraph & g){

        for(unsigned int i = 0; i < path.size(); ++i){
            print_edge(path[i], g);
        }
    }

    void print_paths(EdgeSets paths, const MyGraph & g){

        for(unsigned int i = 0; i < paths.size(); ++i){
            BOOST_LOG_TRIVIAL(info) << "path #" << i;
            print_path(paths[i], g);
        }
    }

    EdgeSet get_edges_from_label(const MyGraph & g, int label){
        EdgeIter ei, ei_end;
        EdgeSet out(g);

        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
            if(g[*ei].label == label)
                out += *ei;
        }
        return out;
    }

    void set_label(EdgeSet p, MyGraph & g, int label){

        EdgeSet::iterator it;
        for (it=p.begin(); it != p.end(); ++it) {
            g[*it].label = label;
        }
    }

    void set_label(Edge e, MyGraph & g, int label){

        g[e].label = label;
    }

    // Assigns a label to all edges
    void set_label_to_all(MyGraph & g, int label){

        EdgeIter ei, ei_end;

        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
            g[*ei].label = label;
        }
    }

    // Assigns a label to edges that are invalid
    void set_label_to_invalid_edges(EdgeSet e_in,
                                    MyGraph & g_in,
                                    MyGraph & g_out,
                                    int label,
                                    bool invert){

        Vertex u;
        Vertex v;

        for(unsigned int i=0; i<e_in.size(); ++i){
            if(invert){
                u = target(e_in[i], g_in);
                v = source(e_in[i], g_in);
            }
            else{
                u = source(e_in[i], g_in);
                v = target(e_in[i], g_in);
            }

            g_out[edge(u,v, g_out).first].label = label;
        }
    }

    // Assigns a label to edges
    void set_label_to_edges(MyGraph & g, EdgeSet es, int label){

        for(unsigned int i=0; i<es.size(); ++i){
            g[es[i]].label = label;
        }
    }

    // Print the whole graph with node/edges ids, weights and labels
    void print_all(const MyGraph & g){

        BOOST_LOG_TRIVIAL(trace) << "------ graph: " << g[graph_bundle].name <<
            " --------";
        BOOST_LOG_TRIVIAL(trace) << "------ edges " << " (n_edges: "
                << num_edges(g) << ") --------";
        EdgeIter ei, ei_end;
        MyGraph::vertex_descriptor v_in;
        MyGraph::vertex_descriptor v_out;

        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
            v_in = source(*ei, g);
            v_out = target(*ei, g);
            BOOST_LOG_TRIVIAL(trace) << "edge: (" << g[v_in].id << "," <<
                g[v_out].id << ") / (" << g[v_in].name << "," << g[v_out].name <<
                ") , id: " << g[*ei].id <<
                ", label: " << g[*ei].label <<
                ", weight: " << g[*ei].weight;
        }

        BOOST_LOG_TRIVIAL(trace) << "------ vertices" << " (n_vertices: "
                << num_vertices(g) << ") --------";
        std::pair<VertexIter, VertexIter> vp;
        for (vp = vertices(g); vp.first != vp.second; ++vp.first)
            BOOST_LOG_TRIVIAL(trace) << "vertex: "
                                    << g[*vp.first].id
                                    << ", name: "
                                    << g[*vp.first].name;
        BOOST_LOG_TRIVIAL(trace) << "--------------";

    }

    // Converts from vertex to edge
    EdgeSet vertpath_to_edgepath(VertexPath path, const MyGraph & g){
        EdgeSet ep(g);
        std::pair<Edge, bool> e_tmp;

        VertexPath::reverse_iterator it;
        for (it=path.rbegin(); it != path.rend()-1; ++it) {
            e_tmp = edge(*it, *std::next(it), g);
            ep += e_tmp.first;
        }
        return ep;
    }

    // To test Bellman-Ford and Dijkstra results
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
                                    int v_id,
                                    const MyGraph & g){

        for(unsigned int i=0; i < p.size(); ++i){
            if(g[source(p[i], g)].id == v_id){
                return i;
            }
        }

        return -1;
    }

    int find_ind_edge_starting_with(EdgeSet p,
                                    Vertex v,
                                    const MyGraph & g){

        for(unsigned int i=0; i < p.size(); ++i){
            if(g[source(p[i], g)].id == g[v].id){
                return i;
            }
        }

        return -1;
    }

    // This function is called when interlacing occurs.
    // One or several edges of p_inter are appended to p.
    EdgeSet append_inter(EdgeSet p,
                         EdgeSet p_inter,
                         Vertex start,
                         Vertex sink,
                         const MyGraph & g){

        // Find label of path p
        int label_p = g[p[0]].label;

        // Find index of start in p_inter
        EdgeSet to_append(g);
        int ind_start = find_ind_edge_starting_with(p_inter, start, g);

        Edge curr_edge = p_inter[ind_start];

        // Append edges of p_inter until:
        // (1) sink is reached
        // (2) target of edge is not part of p (based on label)
        // (3) p_inter is "broken"
        while(true){
            to_append += curr_edge;
            if(g[target(curr_edge, g)].id == g[sink].id)
                break;
            if(g[curr_edge].label == label_p)
                break;
            if(g[p_inter[ind_start]].id_vertex_out !=
               g[p_inter.next(curr_edge)].id_vertex_in)
                break;
            ind_start += 1;
            curr_edge = p_inter[ind_start];
        }

        return p + to_append;
    }

    // Convert an edge set from one graph to another
    // (edge descriptors) are invalidated during remove/add operations
    EdgeSets translate_edge_sets(EdgeSets P,
                                const MyGraph & g_new){

            EdgeSets P_out;
            for(unsigned int i=0; i<P.size(); ++i){
                P[i].convert_to_graph(g_new);
                P_out.push_back(P[i]);
            }
            return P_out;
    }

    Edge translate_edge(Edge e,
                        const MyGraph & g_p,
                        const MyGraph & g){
        std::pair<Edge, bool> e_out;
        e_out = edge(source(e, g_p), target(e, g_p), g);

        if(!e_out.second){
            // cannot be translated! This could be a problem :~(
            //print_edge(e, g_p);
        }
        return e_out.first;
}

    void invert_edge(Edge e,
                     bool inv_label,
                     bool inv_algebraic_sign,
                     MyGraph & g){

        Vertex u = source(e, g);
        Vertex v = target(e, g);

        //print_edge(e, g);

        if(g[e].label > 0){
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
        }
    }

    void invert_edges(EdgeSets edge_sets,
                      bool inv_label,
                      bool inv_algebraic_sign,
                      MyGraph & g){

        for(unsigned int i=0; i<edge_sets.size(); ++i){
            invert_edges(edge_sets[i], inv_label, inv_algebraic_sign, g);
        }
    }

    void invert_edges(EdgeSet edge_path,
                      bool inv_label,
                      bool inv_algebraic_sign,
                      MyGraph & g){

        //print_path(edge_path, g);
        for(unsigned int i=0; i<edge_path.size(); ++i){
            invert_edge(edge_path[i], inv_label, inv_algebraic_sign, g);
        }
    }

    // Returns the total cost of edge sets
    double calc_cost(EdgeSets P, const MyGraph & g){
        double cost;

        for(unsigned int i=0; i<P.size(); ++i){
            for(unsigned int j=0; j<P[i].size(); ++j){
                cost += g[P[i][j]].weight;
            }
        }

        return cost;
    }


    // Convert to python list
    bp::list edgeSets_to_list(EdgeSets P, const MyGraph & g){

        bp::list array;
        for (unsigned int i = 0; i < P.size(); i++){
            boost::python::list temp;
            for (unsigned int j = 0; j < P[i].size(); j++){
                temp.append(g[P[i][j]].id);
            }
            array.append(temp);
        }

        return array;
    }

    EdgeSets augment(EdgeSets P_l,
                     EdgeSet p_inter,
                     Vertex sink_vertex,
                     MyGraph & g,
                     MyGraph & g_c,
                     MyGraph & g_l){

        EdgeSets P_l_plus_1;

        //Loop over P_l. Stop when label = 0, i.e. edge correspond
        // to a "cut"
        // If necessary, look in free edges, i.e. edges in p_inter or leftovers
        // from other interlaced paths
        EdgeSet leftovers(g);
        int label_p = -P_l.size();
        Vertex curr_vertex;
        bool got_cut;
        bool got_non_cont;
        for (auto p : boost::adaptors::reverse(P_l)){
            // Loop over p = P_l[i]
            int ind_edge = 0;
            Edge curr_edge = p[ind_edge];
            int discontinuous_ind;
            while(true){
                got_cut = false;
                got_non_cont = false;
                if(target(curr_edge, g) == sink_vertex)
                    break;

                if(g[curr_edge].label == 0){ // we hit a cut...
                    p.remove_edge(curr_edge);
                    curr_vertex = source(curr_edge, g);
                    p = utils::append_inter(p,
                                            p_inter,
                                            curr_vertex,
                                            sink_vertex,
                                            g);
                    // remove all edges that have been interlaced
                    discontinuous_ind = p.is_discontinuous();
                    while(true){
                        leftovers += p[discontinuous_ind];
                        p.remove_edge(p[discontinuous_ind]);
                        discontinuous_ind = p.is_discontinuous();
                        if(discontinuous_ind == -1)
                            break;
                    }
                    curr_edge = p.back();
                    got_cut = true;
                }
                if(target(curr_edge, g) == sink_vertex)
                    break;
                // no cut but need an edge from p_inter or leftovers
                else if(!p.are_contiguous(curr_edge, p.next(curr_edge)) ||
                        g[target(p.back(), g)].id != g[sink_vertex].id){
                    // unfinished path, look in free edges
                    //utils::print_path(leftovers+p_inter, g);
                    p += build_p_from_self_or_store(curr_edge,
                                                    leftovers + p_inter,
                                                    label_p,
                                                    sink_vertex,
                                                    g);
                    // clear following edges that are not valid anymore
                    discontinuous_ind = p.is_discontinuous();
                    while(true){
                        p.remove_edge(p[discontinuous_ind]);
                        discontinuous_ind = p.is_discontinuous();
                        if(discontinuous_ind == -1)
                            break;
                    }
                    curr_edge = p.back();
                    got_non_cont = true;
                }
                if(target(curr_edge, g) == sink_vertex)
                    break;
                // Nothing special, keep going...
                if(!got_cut || !got_non_cont){
                    curr_edge = p.next(curr_edge);
                }
            }
            P_l_plus_1.push_back(p);
        }

        // Need to invert order of P_l_plus_1
        std::reverse(P_l_plus_1.begin(), P_l_plus_1.end());

        // Building p_ from p_inter and leftovers
        EdgeSet p_ = build_p_from_self_or_store(p_inter[0],
                                                    leftovers+p_inter,
                                                    label_p-1,
                                                    sink_vertex,
                                                    g);

        P_l_plus_1.push_back(p_);

        for(unsigned int i=0; i<P_l_plus_1.size(); ++i)
            utils::set_label_to_edges(g, P_l_plus_1[i], -(i+1));

        return P_l_plus_1;
    }

    EdgeSet build_p_from_self_or_store(Edge start_edge,
                                           EdgeSet store,
                                           int label_p,
                                           Vertex sink_vertex,
                                           const MyGraph & g){

        // build added path p_ from store
        EdgeSet p_(g);
        Edge e;
        MyGraph::out_edge_iterator ei, ei_end;
        p_ += start_edge;
        Vertex curr_vertex = target(start_edge, g);
        int curr_vertex_id = g[curr_vertex].id;
        int ind_store;
        bool got_self = false;
        bool got_store = false;
        //utils::print_path(p_, g);
        while(true){
            ind_store = utils::find_ind_edge_starting_with(store,
                                                            curr_vertex_id,
                                                            g);
            // Check first if we can merge back to myself
            got_self = false;
            got_store = false;

            for (boost::tie(ei, ei_end) = out_edges(curr_vertex, g);
                ei != ei_end; ++ei) {
                e = edge(source(*ei, g),
                            target(*ei, g) , g).first;
                if(g[e].label == label_p){
                    p_ += *ei;
                    got_self = true;
                    curr_vertex = target(p_.back(), g);
                    curr_vertex_id = g[curr_vertex].id;
                    break;
                }
            }
            if(ind_store != -1 && !got_self){ // Add edge from store
                //std::cout << "got from store" << std::endl;
                e = store[ind_store];
                p_ += e;
                got_store = true;
                store = store.remove_edge(e);
                curr_vertex = target(p_.back(), g);
                curr_vertex_id = g[curr_vertex].id;
            }
            if(got_self + got_store == 0){
                //std::cout << "got neither self or store" << std::endl;
                break;
            }
        }
        return p_;
    }
}
