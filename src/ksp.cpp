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
#include "ksp.h"

Ksp::Ksp() {

    new_graph();
}

void Ksp::config(int source_vertex_id,
                 int sink_vertex_id,
                 int a_l_max,
                 std::string source_vertex_name,
                 std::string sink_vertex_name,
                 std::string loglevel,
                 bool a_min_cost){
    /* Define:
       - source_vertex_id: id of source vertex
       - sink_vertex_id: id of source vertex
       - l_max: Max number of paths to compute. -1 for no limit
       - source_vertex_name: Name of source vertex (optional)
       - sink_vertex_name: Name of sink vertex (optional)
       - loglevel: What to print? (see globals.h)
       - min_cost: False for saturation of graph. True to stop at minimum cost
     */

    set_source(source_vertex_id, source_vertex_name);
    set_sink(sink_vertex_id, sink_vertex_name);

    if(a_l_max == -1)
        l_max = std::numeric_limits<int>::max();
    else
        l_max = a_l_max;

    set_loglevel(loglevel);
    min_cost = a_min_cost;

    return_edges = true;
}


Vertex Ksp::add_vertex(int id, std::string str){
    return utils::add_vertex(*G, id, str);
}

void Ksp::remove_edge(int u_id, int v_id){

    Vertex u, v;
    VertexIter wi, wi_end;
    for (boost::tie(wi, wi_end) = vertices(*G); wi != wi_end; ++wi){
        if((*G)[*wi].id == u_id)
            u = *wi;
    }

    for (boost::tie(wi, wi_end) = vertices(*G); wi != wi_end; ++wi){
        if((*G)[*wi].id == v_id)
            v = *wi;
    }

    std::pair<Edge, bool> e = boost::edge(u, v, *G);
    if(e.second)
        boost::remove_edge(u, v, *G);

}

bool Ksp::add_edge(int n0,
                      int n1,
                      double w,
                      int id,
                      std::string str_0,
                      std::string str_1,
                      int label){

    return utils::add_edge(*G, n0, n1, w, id, str_0, str_1, label);
}

void Ksp::set_source(int id, std::string name)
{
    source_vertex = utils::add_vertex(*G, id, name);

}

void Ksp::set_sink(int id, std::string name)
{
    sink_vertex = utils::add_vertex(*G, id, name);
}

int Ksp::num_vertices(){
  return boost::num_vertices(*G);
}

int Ksp::num_edges(){
  return boost::num_edges(*G);
}


void Ksp::new_graph(){

    n_vertices = 0;
    G = new MyGraph(n_vertices);
    (*G)[graph_bundle].name = "defaultName";
    BOOST_LOG_TRIVIAL(debug) << "Creating graph with n_vertices: " <<
        n_vertices;
    BOOST_LOG_TRIVIAL(debug) << "n_edge: " <<
      num_edges();

}

void Ksp::set_loglevel(std::string a_log_level) {

    using namespace boost::log::trivial;
    logSeverityLevel log_level;

    if(a_log_level == "trace")
        log_level = logSeverityLevel::trace;
    else if(a_log_level == "debug")
        log_level = logSeverityLevel::debug;
    else if(a_log_level == "info")
        log_level = logSeverityLevel::info;
    else if(a_log_level == "warning")
        log_level = logSeverityLevel::warning;
    else if(a_log_level == "error")
        log_level = logSeverityLevel::error;
    else
        log_level = logSeverityLevel::info;

    boost::log::core::get()->set_filter
        (
            boost::log::trivial::severity >= log_level
        );
}

bp::list Ksp::run(){
    // Runs the edge-disjoint K-shortest path

    // Make copy of graph for cost transformation
    G_c = new MyGraph(0);

    copy_graph(*G, *G_c);
    (*G_c)[graph_bundle].name = (*G)[graph_bundle].name + "_cost_transform";

    // Store output of shortest-paths
    EdgeSet res_path(*G);
    EdgeSet p_inter(*G_c);
    bool res_ok;
    std::vector<double> res_distance;

    EdgeSets P;
    EdgeSets P_prev;

    BOOST_LOG_TRIVIAL(debug) << "Checking duplicate vertex ids for g";
    BOOST_LOG_TRIVIAL(debug) << utils::has_duplicate_vertex_ids(*G);

    BOOST_LOG_TRIVIAL(debug) << "Checking duplicate edge ids for g";
    BOOST_LOG_TRIVIAL(debug) << utils::has_duplicate_edge_ids(*G);

    BOOST_LOG_TRIVIAL(debug) << "total num edges: "
                            << boost::num_edges(*G);

    BOOST_LOG_TRIVIAL(debug) << "total num nodes: "
                            << boost::num_vertices(*G);

    BOOST_LOG_TRIVIAL(debug) << "num edges leaving source: "
                            << boost::out_degree(source_vertex, *G);

    BOOST_LOG_TRIVIAL(debug) << "num edges entering sink: "
                            << boost::in_degree(sink_vertex, *G);

    double tmp = 0;
    MyGraph::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(*G);
        ei != ei_end; ++ei) {
      tmp += (*G)[*ei].weight;
    }
    BOOST_LOG_TRIVIAL(debug) << "total cost: "
                            << tmp;

    BOOST_LOG_TRIVIAL(debug) << "graph: ";
    utils::print_all(*G);

    BOOST_LOG_TRIVIAL(debug) << "Bellman-Ford...";
    std::tie(res_path, res_ok, res_distance) =
        bellman_ford_shortest_paths(*G);

    if(!res_ok){
        BOOST_LOG_TRIVIAL(error) << "Couldn't compute a single path!";
    }
    else{
        BOOST_LOG_TRIVIAL(debug) << "...ok";
        P.push_back(res_path);
        cost = utils::calc_cost(P, *G);
        BOOST_LOG_TRIVIAL(info) << "l: " << 0
                                << ", cost: " << cost;

        utils::print_paths(P, *G);
    }

    utils::invert_edges(P,
                        true,
                        true,
                        *G_c);

    BOOST_LOG_TRIVIAL(debug) << "inverted edges on g_c";
    P_prev = P;

    for(int l = 1; l < l_max; ++l){

        EdgeSet p_cut(*G);

        cost_transform(res_distance, *G_c);
        BOOST_LOG_TRIVIAL(debug) << "done cost_transform";
        utils::print_all(*G_c);
        BOOST_LOG_TRIVIAL(debug) << "Dijkstra...";
        std::tie(p_inter, res_ok, res_distance) =
            dijkstra_shortest_paths(*G_c, source_vertex);

        // p_cut are invalid edges, i.e. will be excluded from augmentation set
        p_cut = p_inter.convert_to_graph(*G);
        BOOST_LOG_TRIVIAL(debug) << "p_inter:";
        utils::print_path(p_inter, *G);

        // convert p_cut to valid edges on *G
        Edge e;
        Vertex u, v;
        EdgeSet p_cut_g(*G);
        for(unsigned int i=0; i<p_cut.size(); ++i){
            e = p_cut[i];
            u = source(e, *G_c);
            v = target(e, *G_c);
            p_cut_g += edge(v, u ,*G).first;
        }

        if(res_ok){
            BOOST_LOG_TRIVIAL(debug) << "ok...";

            // If interlacing path didn't cut edges of previous solution,
            // add path as is, else augment
            if(p_cut.size() > 0){
                utils::set_label_to_invalid_edges(p_cut, *G_c, *G, 0, true);
                BOOST_LOG_TRIVIAL(debug) << "num invalid edges (cut) : "
                                         << utils::num_edges_with_label(*G, 0);
                BOOST_LOG_TRIVIAL(debug) << "p_cut: ";
                utils::print_path(p_cut, *G_c);
                utils::print_all(*G);

                BOOST_LOG_TRIVIAL(debug) << "Augmenting";
                P = utils::augment(P,
                                   p_cut_g,
                                   p_inter,
                                   sink_vertex,
                                   *G);
            }
            else{
                P.push_back(p_inter);
                utils::set_label(p_inter, *G, -P.size());
            }
            new_cost = utils::calc_cost(P, *G);
            BOOST_LOG_TRIVIAL(info) << "l: " << l
                                    << ", cost: " << new_cost;

            utils::print_paths(P, *G);
            bool dup_edges;
            dup_edges = utils::has_duplicate_edge_ids(utils::flatten(P, *G), *G);
            BOOST_LOG_TRIVIAL(debug) << "duplicate edges in solution: "
                                     << dup_edges;
            utils::invert_edges(utils::translate_edge_sets(P,*G_c),
                                true,
                                true,
                                *G_c);

            // set label=1 and invert back edges of p_cut
            Edge e;
            Vertex u, v;
            for(unsigned int i=0; i<p_cut.size(); ++i){
                e = p_cut[i];
                u = source(e, *G_c);
                v = target(e, *G_c);
                utils::set_label(p_cut, *G_c, 1);
                utils::invert_edge(edge(u,v,*G_c).first,
                                   false,
                                   false,
                                   *G_c);
            }
            utils::set_label_to_invalid_edges(p_cut, *G_c, *G, 1, true);

            utils::print_all(*G_c);
            utils::print_all(*G);

            if(min_cost && (new_cost > cost)){
                BOOST_LOG_TRIVIAL(info) << "Reached minimum at l= " << l-1;
                return utils::edgeSets_to_edges_list(P_prev,
                                                     *G);

            }
            else{
              P_prev = P;
              cost = new_cost;

            }
        }
        else{
            BOOST_LOG_TRIVIAL(info) << "Stopped at l= " << l-1;

            break;
        }
    }

    // Re-initialize edge labels
    BOOST_LOG_TRIVIAL(debug) << "setting all labels to 1";
    utils::set_label_to_all(*G, 1);

    return utils::edgeSets_to_edges_list(P,
                                        *G);
}

std::tuple<EdgeSet, bool, std::vector<double>>
Ksp::dijkstra_shortest_paths(const MyGraph & g, Vertex source_vertex){

    EdgeSet out_path(g);

    // init the distance
    std::vector<double> distances( num_vertices());

    std::vector<std::size_t> predecessors(num_vertices());

    // call to the algorithm
    //bool r = boost::dijkstra_shortest_paths(
    boost::dijkstra_shortest_paths(
        g,
        source_vertex,
        weight_map(get(&EdgeProperty::weight, g)).
        distance_map(make_iterator_property_map(distances.begin(),
                                                       get(vertex_index,g))).
        predecessor_map(boost::make_iterator_property_map(predecessors.begin(),
                                                          get(vertex_index,g)))
        );

    //print_distances(distance, *G);

    VertexPath shortest;
    bool r;

    std::tie(shortest, r) = utils::pred_to_path(predecessors, g,
                                                source_vertex, sink_vertex);

    if(!r){
        BOOST_LOG_TRIVIAL(info) << "Couldn't reach sink node";
        out_path = EdgeSet(g);
        r = false;
    }
    else{
        out_path = utils::vertpath_to_edgepath(shortest, g);
        //utils::print_path(out_path, g);
        r = true;
    }

    return std::make_tuple(out_path, r, distances);
}

std::tuple<EdgeSet, bool, std::vector<double>>
Ksp::bellman_ford_shortest_paths(const MyGraph & g){

    EdgeSet out_path(g);
    std::pair<VertexIter, VertexIter> vp;
    auto v_index = get(boost::vertex_index, g);
    auto weight = boost::make_transform_value_property_map(
        std::mem_fn(&EdgeProperty::weight),
        get(boost::edge_bundle, g));

    // init
    std::vector<Vertex> predecessors(num_vertices());
    //for (vp = vertices(*G); vp.first != vp.second; ++vp.first)
    //    predecessors[v_index[*vp.first]] = v_index[*vp.first];

    std::vector<double> distances(num_vertices(),
                                  (std::numeric_limits<double>::max)());
    distances[v_index[source_vertex]] = 0;

    bool r = boost::bellman_ford_shortest_paths(
        *G,
        num_vertices(),
        weight_map(weight).
        distance_map(make_iterator_property_map(distances.begin(),
                                                v_index)).
        predecessor_map(make_iterator_property_map(predecessors.begin(),
                                                   v_index))
        );

    VertexPath shortest;

    if (r){
        std::tie(shortest, std::ignore) = utils::pred_to_path(predecessors,
                                                              g,
                                                              source_vertex,
                                                              sink_vertex);
        out_path = utils::vertpath_to_edgepath(shortest, g);
    }
    else{
        BOOST_LOG_TRIVIAL(info) << "negative cycle";
        out_path = EdgeSet(g);
    }

    return std::make_tuple(out_path, r, distances);
}

void Ksp::cost_transform(const std::vector<double> & distance,
                         MyGraph & g_out){

    EdgeIter ei, ei_end;
    double s_i;
    double s_j;
    double w_ij;
    for (boost::tie(ei, ei_end) = edges(g_out); ei != ei_end; ++ei){
        s_i = distance[source(*ei, g_out)];
        s_j = distance[target(*ei, g_out)];
        w_ij = g_out[*ei].weight;
        g_out[*ei].weight = w_ij + s_i - s_j;
    }

}

void Ksp::set_label_all_edges(int label){

    utils::set_label_to_all(*G, label);

}
