#include "boost_ksp.h"

Vertex ksp::add_vertex(int id, std::string str){
    return utils::add_vertex(*G, id, str);
}

bool ksp::add_edge(int n0,
                      int n1,
                      double w,
                      int id,
                      std::string str_0,
                      std::string str_1,
                      int label){

    return utils::add_edge(*G, n0, n1, w, id, str_0, str_1, label);
}

void ksp::set_source(int id, std::string str)
{
    source_vertex = utils::add_vertex(*G, id, str);

}

void ksp::set_sink(int id, std::string str)
{
    sink_vertex = utils::add_vertex(*G, id, str);
}

void ksp::new_graph(int a_n_vertices=0){

    n_vertices = a_n_vertices;
    std::cout << "Creating graph with n_vertices: " <<
              n_vertices << std::endl;
    G = new MyGraph(n_vertices);
    (*G)[graph_bundle].name = "defaultName";

}

shared_ptr<ksp> ksp::create() {

    return shared_ptr<ksp>(new ksp);
}

bool ksp::do_ksp(){

    //Copy G to G_c and work with latter

    G_c = new MyGraph(0);
    G_l = new MyGraph(0);

    copy_graph(*G, *G_c);
    (*G_c)[graph_bundle].name = (*G)[graph_bundle].name + "_cost_transform";

    copy_graph(*G, *G_l);
    (*G_l)[graph_bundle].name = (*G)[graph_bundle].name + "_inv_edges";
    //copy_graph(*G, *G_l);

    // Store output of shortest-paths
    ShortestPathRes res;
    EdgeSet res_path;
    bool res_ok;
    std::vector<double> res_distance;

    EdgeSets P;

    //int l_max = std::numeric_limits<int>::max();
    int l_max = 3;

    std::cout << "Bellman-Ford" << std::endl;
    res = bellman_ford_shortest_paths(*G);
    std::tie(res_path, res_ok, res_distance) = res;
    utils::print_path(res_path, *G);

    P.push_back(res_path);

    for(int l = 1; l < l_max; ++l){

        std::cout << "l: " << l << std::endl;

        if(l != 0){
            // Check costs for minima
        }

        for(unsigned int i=0; i<P.size(); ++i)
            utils::set_label_on_path(P[i], -(i+1), *G);

        //std::cout << "before invert edges" << std::endl;
        //utils::print_all(*G_l);
        //std::cout << "after invert edges" << std::endl;
        utils::invert_edges(P, true, true, *G, *G_l);
        utils::invert_edges(P, true, true, *G, *G_c);
        //utils::print_all(*G_l);
        //utils::print_all(*G_c);
        cost_transform(res_distance, *G_l, *G_c);
        std::cout << "done cost_transform" << std::endl;
        utils::print_all(*G_c);
        std::cout << "Dijkstra" << std::endl;
        res = dijkstra_shortest_paths(*G_c);
        std::tie(res_path, res_ok, res_distance) = res;
        utils::print_path(res_path, *G);
        std::cout << "Augment" << std::endl;

        P = augment(P,
                    res_path,
                    *G,
                    *G_c,
                    *G_l);
        std::cout << "Solution at l= "<< l << std::endl;
        utils::print_paths(P, *G);
        utils::print_all(*G);

    }
    return true;
}

ShortestPathRes ksp::dijkstra_shortest_paths(const MyGraph & g){

    EdgeSet out_path;
    ShortestPathRes out;

    // init the distance
    std::vector<double> distances( num_vertices(g));

    std::vector<std::size_t> predecessors(num_vertices(g));

    // call to the algorithm
    //bool r = boost::dijkstra_shortest_paths(
    bool r = true;
    boost::dijkstra_shortest_paths(
        g,
        source_vertex,
        weight_map(get(&EdgeProperty::weight, g)).
        distance_map(make_iterator_property_map(distances.begin(),
                                                       get(vertex_index,g))).
        predecessor_map(boost::make_iterator_property_map(predecessors.begin(),
                                                          get(vertex_index,g)))
        );


    if (r){
        //print_distances(distance, *G);
        VertexPath shortest = utils::pred_to_path(predecessors, g,
                                           source_vertex, sink_vertex);
        //print_path(shortest, *G);
        out_path = utils::vertpath_to_edgepath(shortest, g);

    }
    else{
        std::cout << "negative cycle" << std::endl;
        out_path = EdgeSet();
    }

    out = std::make_tuple(out_path, r, distances);

    return out;
}

ShortestPathRes  ksp::bellman_ford_shortest_paths(const MyGraph & g){
    ShortestPathRes   out;
    EdgeSet out_path;
    std::pair<VertexIter, VertexIter> vp;
    auto v_index = get(boost::vertex_index, g);
    auto weight = boost::make_transform_value_property_map(
        std::mem_fn(&EdgeProperty::weight),
        get(boost::edge_bundle, g));

    // init
    std::vector<Vertex> predecessors(num_vertices(g));
    //for (vp = vertices(*G); vp.first != vp.second; ++vp.first)
    //    predecessors[v_index[*vp.first]] = v_index[*vp.first];

    std::vector<double> distances(num_vertices(g),
                                  (std::numeric_limits<double>::max)());
    distances[v_index[source_vertex]] = 0;

    bool r = boost::bellman_ford_shortest_paths(
        *G,
        num_vertices(g),
        weight_map(weight).
        distance_map(make_iterator_property_map(distances.begin(),
                                                v_index)).
        predecessor_map(make_iterator_property_map(predecessors.begin(),
                                                   v_index))
        );

    if (r){
        //print_dist_pred(distances, predecessors, *G);

        VertexPath shortest = utils::pred_to_path(predecessors, g,
                                           source_vertex, sink_vertex);
        out_path = utils::vertpath_to_edgepath(shortest, g);

    }
    else{
        std::cout << "negative cycle" << std::endl;
        out_path = EdgeSet();
    }

    out = std::make_tuple(out_path, r, distances);

    return out;
}

void ksp::cost_transform(const std::vector<double> & distance,
                         const MyGraph & g_in,
                         MyGraph & g_out){

    EdgeIter ei, ei_end;
    double s_i;
    double s_j;
    double w_ij;
    for (boost::tie(ei, ei_end) = edges(g_out); ei != ei_end; ++ei){
        s_i = distance[source(*ei, g_in)];
        s_j = distance[target(*ei, g_in)];
        w_ij = g_out[*ei].weight;
        g_out[*ei].weight = w_ij + s_i - s_j;
    }
}

EdgeSets ksp::augment(EdgeSets P_l,
                      EdgeSet p_inter,
                      MyGraph & g,
                      MyGraph & g_c,
                      MyGraph & g_l){
    // p_inter is defined on g_c
    // P_l is defined on g

    EdgeSets P_l_plus_1;
    EdgeSet p_cut; // edges of label -1 occupied by p_inter

    Vertex v_in;
    Vertex v_out;

    //Populate p_cut.
    for(unsigned int i=0; i < p_inter.size(); ++i){
        v_in = source(p_inter[i], g);
        v_out = target(p_inter[i], g);
        if(!edge(v_in, v_out, g).second){
            std::cout << "cut at edge: ("  <<
                g[v_in].name << "," <<
                g[v_out].name << ")" << std::endl;
            p_cut.push_back(p_inter[i]);
        }
    }

    // Erase edges of p_cut from p_inter
    for(unsigned int i=0; i < p_cut.size(); ++i)
        p_inter = utils::remove_edge_from_set(p_cut[i],
                                        p_inter,
                                       g,
                                       false);
    std::cout << "removed edges of p_cut from p_inter" << std::endl;

    //Loop over P_l and remove edges belonging to p_cut
    EdgeSet leftovers;
    std::tuple<EdgeSet, EdgeSet, EdgeSet, Edge> res_append;
    int label_p = -P_l.size();
    for (auto p : boost::adaptors::reverse(P_l)){
        // Loop over p = P_l[i]
        int ind_edge = 0;
        Edge curr_edge = p[ind_edge];
        while(target(curr_edge, g) != sink_vertex){
            utils::print_edge(curr_edge, g);
            //print_path(p_cut, *G);
            if(utils::edge_is_in_set(curr_edge, p_cut, g, g_c, true)){
                std::cout << "edge is in p_cut" << std::endl;
                p = utils::remove_edge_from_set(curr_edge,
                                         p,
                                         g,
                                         false);
                res_append = utils::append_inter(p,
                                                 p_inter,
                                                 p_cut,
                                                 leftovers,
                                                 source(curr_edge, g),
                                                 sink_vertex,
                                                 g,
                                                 g_c);
                std::tie(p, p_inter, leftovers, curr_edge) = res_append;
                std::cout << "p" << std::endl;
                utils::print_path(p, g);
                std::cout << "p_inter" << std::endl;
                utils::print_path(p_inter, g);
                std::cout << "leftovers" << std::endl;
                utils::print_path(leftovers, g);
            }
            else{
                ind_edge += 1;
                curr_edge = p[ind_edge];
            }
        }
        P_l_plus_1.push_back(p);
    }

    // p_cut are re-established on g_c
    utils::invert_edges(p_cut, true, true, g_c, g_c);
    utils::invert_edges(p_cut, true, true, g_l, g_l);

    // Need to invert order of P_l_plus_1
    std::reverse(P_l_plus_1.begin(), P_l_plus_1.end());

    // build added path p_ from p_inter and leftovers
    EdgeSet p_;
    MyGraph::out_edge_iterator ei, ei_end;
    Edge e = p_inter[0];
    e = edge(source(e, g_c), target(e, g_c) , g).first;
    p_.push_back(e);
    Vertex curr_vertex = target(p_[0], g);
    int ind;
    while(true){
        if(g[curr_vertex].id == g[sink_vertex].id)
            break;
        // Search in p_inter for next edge
        ind = utils::find_ind_edge_starting_with(p_inter,
                                                 curr_vertex,
                                                 g);
        if(ind == -1){ // Search in leftovers
            for (boost::tie(ei, ei_end) = out_edges(curr_vertex, g);
                ei != ei_end; ++ei) {
                e = edge(source(*ei, g_c),
                         target(*ei, g_c) , g).first;
                 if(g[e].label == label_p){
                    p_.push_back(*ei);
                    break;
                }
            }
        }
        else
            p_.push_back(p_inter[ind]);
        //std::cout << "curr_vertex name: " << (*G)[curr_vertex].name << std::endl;
        curr_vertex = target(p_.back(), g);
    }
    p_ = utils::translate_edge_set(p_, g_c, g);
    P_l_plus_1.push_back(p_);

    //std::cout << "setting label=1 on edges:" << std::endl;
    //utils::set_label_to_all(g, 1);
    utils::print_all(g);
    for(unsigned int i=0; i<P_l_plus_1.size(); ++i)
        utils::set_label_to_edges(g, P_l_plus_1[i], -(i+1));
    utils::print_all(g);


    //utils::print_all(g_l);
    //utils::print_all(g_c);
    //std::cout << "inverting edges on g_l and g_c:" << std::endl;
    //utils::invert_edges(p_cut, false, false, g_c, g_c);
    //utils::invert_edges(p_cut, false, false, g_l, g_l);
    //utils::print_all(g_l);
    //utils::print_all(g_c);


    return P_l_plus_1;
}
