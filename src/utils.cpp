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

        std::pair<Edge, bool> e = boost::add_edge(v0, v1, g);

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

    void print_edge(Edge e, const MyGraph & g){
        BOOST_LOG_TRIVIAL(debug) << "(" << g[source(e, g)].name
                                 << ","
                                 << g[target(e, g)].name << ")"
                                 << "/"
                                 << "(" << g[source(e, g)].id
                                 << "," << g[target(e, g)].id << ") "
                                 << "label: " << g[e].label;
    }

    void print_path(EdgeSet path, const MyGraph & g){

        for(unsigned int i = 0; i < path.size(); ++i){
            print_edge(path[i], g);
        }
    }

    void print_paths(EdgeSets paths, const MyGraph & g){

        for(unsigned int i = 0; i < paths.size(); ++i){
            BOOST_LOG_TRIVIAL(debug) << "path #" << i;
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

            if(!edge(u,v, g_out).second)
                BOOST_LOG_TRIVIAL(debug) <<
                  "set_label_to_invalid_edges: Edge doesnt exist!!!";

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

      BOOST_LOG_TRIVIAL(debug) << "find_ind_edge_starting_with (int): " << v_id;
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
    // One edge of p_app is returned
    Edge append_edge(EdgeSet p_app,
                    Vertex start,
                    const MyGraph & g){

        // Find index of start in p_app
        unsigned int ind_start = find_ind_edge_starting_with(p_app,
                                                             g[start].id,
                                                             g);

        return p_app[ind_start];
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
        double cost = 0;

        for(unsigned int i=0; i<P.size(); ++i){
            for(unsigned int j=0; j<P[i].size(); ++j){
                cost += g[P[i][j]].weight;
            }
        }

        return cost;
    }

    // Convert to python list with edges ids
    bp::list edgeSets_to_edges_list(EdgeSets P, const MyGraph & g){

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
                     EdgeSet p_cut,
                     EdgeSet p_inter,
                     Vertex sink_vertex,
                     MyGraph & g){

        // initialize P_l_plus_1 with empty sets
        EdgeSets P_l_plus_1;
        EdgeSet P_l_plus_1_flat(g);

        // This stores edges of last run to pick from
        EdgeSet P_l_clean = flatten(P_l, g) - p_cut;
        p_inter -= p_cut;

        // Building p_ from p_inter and store
        BOOST_LOG_TRIVIAL(debug) << "Building p_ from p_inter and last solution";
        EdgeSet p_(g);

        Edge curr_edge = p_inter[0];
        Vertex curr_vertex;
        p_ += curr_edge;
        print_edge(curr_edge, g);
        while(target(curr_edge, g) != sink_vertex){
            curr_vertex = target(curr_edge, g);
            if(p_inter.has_out_vertex(curr_vertex)){
                BOOST_LOG_TRIVIAL(debug) << "append inter";
                curr_edge = append_edge(p_inter,
                                        curr_vertex,
                                        g);
                // remove all edges that have been interlaced
                BOOST_LOG_TRIVIAL(debug) << "remove interlaced edges";

                // clean p_inter
                p_inter -= curr_edge;

            }
            else{
                BOOST_LOG_TRIVIAL(debug) << "pick from last solution set";
                curr_edge = last_out_edge(curr_vertex, P_l_clean, g);
                P_l_clean -= curr_edge;

            }
            p_ += curr_edge;
            curr_edge = p_.back();
            print_edge(curr_edge, g);
        }

        P_l_clean.sort_descend_labels();

        set_label(p_, g, -(P_l.size()-1));

        int label_p;
        EdgeSet edges_self(g);
        for(auto p : boost::adaptors::reverse(P_l)){
            label_p = g[p[0]].label; // we are tracking this label
            Edge curr_edge = p[0];
            EdgeSet p_new(g);
            p_new += curr_edge;
            P_l_clean -= curr_edge;
            while(target(curr_edge, g) != sink_vertex){
                BOOST_LOG_TRIVIAL(debug) << "are P_l_clean labels sorted: "
                                            << P_l_clean.are_label_sorted();
                BOOST_LOG_TRIVIAL(debug) << "label: " << label_p;
                curr_vertex = target(curr_edge, g);
                edges_self = out_edges_with_label(curr_vertex,
                                                  label_p,
                                                  P_l_clean,
                                                  g);
                //append inter or branch to other path
                if(edges_self.size() == 0){
                    if(p_inter.has_out_vertex(curr_vertex)){
                        BOOST_LOG_TRIVIAL(debug) << "append inter";
                        curr_edge = append_edge(p_inter,
                                                curr_vertex,
                                                g);
                        // remove all edges that have been interlaced
                        BOOST_LOG_TRIVIAL(debug) << "remove interlaced edges";

                        // clean p_inter
                        p_new += curr_edge;
                        p_inter -= curr_edge;
                        P_l_clean -= curr_edge;
                    }
                    else{
                        BOOST_LOG_TRIVIAL(debug) << "branch to another path";
                        curr_edge = last_out_edge(curr_vertex, P_l_clean, g);
                        p_new += curr_edge;
                        P_l_clean -= curr_edge;
                        BOOST_LOG_TRIVIAL(debug) << "end branch to another path";
                    }
                    curr_edge = p_new.back();
                    //set_label(curr_edge, g, label_p);
                }
                else{ //follow myself
                    curr_edge = edges_self[0];
                    p_new += curr_edge;
                    P_l_clean -= curr_edge;
                }

                print_edge(curr_edge, g);
            }
            set_label(p_new, g, label_p);
            P_l_plus_1.insert(P_l_plus_1.begin(), p_new);
        }

        BOOST_LOG_TRIVIAL(debug) << "done augmenting l paths";
        P_l_plus_1.push_back(p_);

        return P_l_plus_1;
    }

    bool has_duplicate_vertex_ids(const MyGraph & g){

        std::vector<int> ids;
        VertexIter vi, vi_end;

        for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi){
            ids.push_back(g[*vi].id);
        }

        auto it = std::unique( ids.begin(), ids.end() );
        return !(it == ids.end() );
    }

    bool has_duplicate_edge_ids(EdgeSet p, const MyGraph & g){

        std::vector<int> ids;
        EdgeSet::iterator ei;

        for(ei=p.begin(); ei!=p.end(); ++ei){
            ids.push_back(g[*ei].id);
        }

        auto it = std::unique( ids.begin(), ids.end() );
        return !(it == ids.end() );
    }

    bool has_duplicate_edge_ids(const MyGraph & g){

        std::vector<int> ids;
        EdgeIter ei, ei_end;

        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
            ids.push_back(g[*ei].id);
        }

        auto it = std::unique( ids.begin(), ids.end() );
        return !(it == ids.end() );
    }

    int num_edges_with_label(const MyGraph & g, int label){

            std::vector<int> ids;
            EdgeIter ei, ei_end;

            for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
            if(g[*ei].label == label)
                ids.push_back(g[*ei].id);
            }

            return ids.size();
        }

    EdgeSet out_edges_with_label(Vertex v, int label, const MyGraph & g){
        // return out_edges on graph g at vertex v with label label.

        EdgeSet out(g);
        MyGraph::out_edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = out_edges(v, g);
             ei != ei_end; ++ei) {
            if(g[*ei].label == label)
                out += *ei;

        }

        return out;
    }

    EdgeSet out_edges_with_label(Vertex v, int label, EdgeSet & s,
                                 const MyGraph & g){
        // return out_edges on graph g at vertex v with label label.

        EdgeSet out(g);
        EdgeSet::iterator it;
        for(it=s.begin(); it!= s.end(); ++it){
            if((g[*it].label == label) && (source(*it, g) == v))
                out += *it;
        }

        return out;
    }

    EdgeSet out_edges_with_neg_label(Vertex v, const MyGraph & g){
        // return out_edges on graph g at vertex v with negative label.

        EdgeSet out(g);
        MyGraph::out_edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = out_edges(v, g);
             ei != ei_end; ++ei) {
            if(g[*ei].label < 0){
                out += *ei;
            }

        }

        return out;
    }

    Edge first_out_edge(Vertex v, EdgeSet& s,
                                           const MyGraph & g){
        // return edges on graph g at vertex v. Returns first instance.
        EdgeSet::iterator it;
        for(it=s.begin(); it!= s.end(); ++it){
            if(g[source(*it, g)].id == g[v].id){
                return *it;
            }

        }

        BOOST_LOG_TRIVIAL(debug) << "first_out_edge: found no edge with in vertex: " << g[v].id;
        return *it;
    }

    Edge last_out_edge(Vertex v, EdgeSet& s,
                                           const MyGraph & g){
        // return edges on graph g at vertex v. Returns last instance.
        EdgeSet::reverse_iterator it;
        for(it=s.rbegin(); it!= s.rend(); ++it){
            if(g[source(*it, g)].id == g[v].id){
                BOOST_LOG_TRIVIAL(debug) << "last_out_edge: return label: " << g[*it].label;
                return *it;
            }

        }
        BOOST_LOG_TRIVIAL(debug) << "last_out_edge: found no edge with in vertex: " << g[v].id;

        return *it;
    }

    EdgeSet out_edges_with_neg_label(Vertex v, EdgeSet& s, const MyGraph & g){
        // return edges of set s on graph g at vertex v with negative label.

        EdgeSet out(g);
        EdgeSet::iterator it;
        for(it=s.begin(); it!= s.end(); ++it){
            if((g[*it].label < 0) && (g[source(*it, g)].id == g[v].id)){
                out += *it;
            }

        }

        return out;
    }

    EdgeSet flatten(EdgeSets & P, const MyGraph& g){

        EdgeSet out = EdgeSet(g);
        for(unsigned int i=0; i<P.size(); ++i){
            for(unsigned int j=0; j<P[i].size(); ++j)
                out += P[i][j];
        }
        return out;
    }
}
