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

}

shared_ptr<ksp> ksp::create() {

    return shared_ptr<ksp>(new ksp);
}

bool ksp::do_ksp(){

    // Store output of shortest-paths
    ShortestPathRes res;
    EdgeSet res_path;
    bool res_ok;
    std::vector<double> res_distance;

    EdgeSets P;

    //int l_max = std::numeric_limits<int>::max();
    int l_max = 3;

    std::cout << "Bellman-Ford" << std::endl;
    res = bellman_ford_shortest_paths();
    std::tie(res_path, res_ok, res_distance) = res;
    utils::print_path(res_path, *G);

    P.push_back(res_path);

    for(int l = 1; l < l_max; ++l){

        std::cout << "l: " << l << std::endl;

        if(l != 0){
            // Check costs for minima
        }

        utils::invert_edges(P, true, true, *G);
        cost_transform(res_distance);
        std::cout << "done cost_transform" << std::endl;
        utils::print_all(*G);
        std::cout << "Dijkstra" << std::endl;
        res = dijkstra_shortest_paths();
        std::tie(res_path, res_ok, res_distance) = res;
        utils::print_path(res_path, *G);
        std::cout << "Augment" << std::endl;

        P = augment(P, res_path);
        std::cout << "Solution at l= "<< l << std::endl;
        utils::print_paths(P, *G);
        utils::print_all(*G);

    }
    return true;
}

ShortestPathRes ksp::dijkstra_shortest_paths(){

    EdgeSet out_path;
    ShortestPathRes out;

    // init the distance
    std::vector<double> distances( num_vertices(*G));

    std::vector<std::size_t> predecessors(num_vertices(*G));

    // call to the algorithm
    //bool r = boost::dijkstra_shortest_paths(
    bool r = true;
    boost::dijkstra_shortest_paths(
        *G,
        source_vertex,
        weight_map(get(&MyEdge::weight, *G)).
        distance_map(make_iterator_property_map(distances.begin(),
                                                       get(vertex_index,*G))).
        predecessor_map(boost::make_iterator_property_map(predecessors.begin(),
                                                          get(vertex_index,*G)))
        );


    if (r){
        //print_distances(distance, *G);
        VertexPath shortest = utils::pred_to_path(predecessors, *G,
                                           source_vertex, sink_vertex);
        //print_path(shortest, *G);
        out_path = utils::vertpath_to_edgepath(shortest, *G);

    }
    else{
        std::cout << "negative cycle" << std::endl;
        out_path = EdgeSet();
    }

    out = std::make_tuple(out_path, r, distances);

    return out;
}

ShortestPathRes  ksp::bellman_ford_shortest_paths(){
    ShortestPathRes   out;
    EdgeSet out_path;
    std::pair<VertexIter, VertexIter> vp;
    auto v_index = get(boost::vertex_index, *G);
    auto weight = boost::make_transform_value_property_map(
        std::mem_fn(&MyEdge::weight),
        get(boost::edge_bundle, *G));

    // init
    std::vector<Vertex> predecessors(num_vertices(*G));
    //for (vp = vertices(*G); vp.first != vp.second; ++vp.first)
    //    predecessors[v_index[*vp.first]] = v_index[*vp.first];

    std::vector<double> distances(num_vertices(*G),
                                  (std::numeric_limits<double>::max)());
    distances[v_index[source_vertex]] = 0;

    bool r = boost::bellman_ford_shortest_paths(
        *G,
        num_vertices(*G),
        weight_map(weight).
        distance_map(make_iterator_property_map(distances.begin(),
                                                v_index)).
        predecessor_map(make_iterator_property_map(predecessors.begin(),
                                                   v_index))
        );

    if (r){
        //print_dist_pred(distances, predecessors, *G);

        VertexPath shortest = utils::pred_to_path(predecessors, *G,
                                           source_vertex, sink_vertex);
        out_path = utils::vertpath_to_edgepath(shortest, *G);

    }
    else{
        std::cout << "negative cycle" << std::endl;
        out_path = EdgeSet();
    }

    out = std::make_tuple(out_path, r, distances);

    return out;
}

void ksp::cost_transform(const std::vector<double> & distance){

    EdgeIter ei, ei_end;
    double s_i;
    double s_j;
    double w_ij;
    for (boost::tie(ei, ei_end) = edges(*G); ei != ei_end; ++ei){
        s_i = distance[source(*ei, *G)];
        s_j = distance[target(*ei, *G)];
        w_ij = (*G)[*ei].weight;
        (*G)[*ei].weight = w_ij + s_i - s_j;
    }
}

EdgeSets ksp::augment(EdgeSets P_l, EdgeSet p_inter){
    EdgeSets P_l_plus_1;
    EdgeSet p_cut; // edges of label -1 occupied by p_inter

    Vertex v_in;
    Vertex v_out;

    //Populate p_cut
    for(unsigned int i=0; i < p_inter.size(); ++i){
        v_in = source(p_inter[i], *G);
        v_out = target(p_inter[i], *G);
        if((*G)[p_inter[i]].label == -1){
            std::cout << "cut at edge: ("  <<
                (*G)[v_in].name << "," <<
                (*G)[v_out].name << ")" << std::endl;
            p_cut.push_back(p_inter[i]);
        }
    }

    // Erase edges of p_cut from p_inter
    for(unsigned int i=0; i < p_cut.size(); ++i)
        p_inter = utils::remove_edge_from_set(p_cut[i],
                                        p_inter,
                                       *G,
                                       false);
    std::cout << "removed edges of p_cut from p_inter" << std::endl;


    //Loop over P_l and remove edges belonging to p_cut
    EdgeSet leftovers;
    std::tuple<EdgeSet, EdgeSet, EdgeSet, Edge> res_append;
    for(unsigned int i=0; i < P_l.size(); ++i){
        EdgeSet p = P_l[i];

        // Loop over p = P_l[i]
        int ind_edge = 0;
        Edge curr_edge = p[ind_edge];
        while(target(curr_edge, *G) != sink_vertex){
            //print_edge(curr_edge, *G);
            //print_path(p_cut, *G);
            if(utils::edge_is_in_set(curr_edge, p_cut, *G, true)){
                //std::cout << "edge is in p_cut" << std::endl;
                p = utils::remove_edge_from_set(curr_edge,
                                         p,
                                         *G,
                                         false);
                res_append = utils::append_inter(p,
                                                 p_inter,
                                                 p_cut,
                                                 leftovers,
                                                 (*G)[source(curr_edge, *G)],
                                                 (*G)[sink_vertex],
                                                 *G);

                std::tie(p, p_inter, leftovers, curr_edge) = res_append;
                //std::cout << "p" << std::endl;
                //print_path(p, *G);
                //std::cout << "p_inter" << std::endl;
                //print_path(p_inter, *G);
                //std::cout << "leftovers" << std::endl;
                //print_path(leftovers, *G);
            }
            else{
                ind_edge += 1;
                curr_edge = p[ind_edge];
            }

        }

        P_l_plus_1.push_back(p);
    }

    // invert directions (and algebraic signs?) of all paths in P_l
    //print_all(*G);
    //std::cout << "getting edges with label = -1" << std::endl;
    EdgeSet to_invert = utils::get_edges_from_label(*G, -1);
    //std::cout << "inverting paths in augment" << std::endl;
    for(unsigned int i = 0; i<to_invert.size(); ++i){
        utils::invert_edge(to_invert[i],
                           false,
                           true,
                           *G);
    }

    // build added path p_ from p_inter and leftovers
    EdgeSet p_;
    p_.push_back(p_inter[0]);
    Vertex curr_vertex = target(p_[0], *G);
    MyGraph::out_edge_iterator ei, ei_end;
    Edge e;
    int ind;
    while(true){
        if((*G)[curr_vertex].id == (*G)[sink_vertex].id)
            break;
        // Search in p_inter for next edge
        ind = utils::find_ind_edge_starting_with(p_inter,
                                                 (*G)[curr_vertex],
                                                 (*G));
        if(ind == -1){ // Search in leftovers
            for (boost::tie(ei, ei_end) = out_edges(curr_vertex, *G);
                ei != ei_end; ++ei) {
                 if((*G)[*ei].label == -1){
                    p_.push_back(*ei);
                    break;
                }
            }
        }
        else
            p_.push_back(p_inter[ind]);

        //std::cout << "curr_vertex name: " << (*G)[curr_vertex].name << std::endl;
        curr_vertex = target(p_.back(), *G);
    }
    P_l_plus_1.push_back(p_);

    //std::cout << "setting label=1 on edges:" << std::endl;
    utils::set_label_to_all(*G, 1);
    utils::print_all(*G);

    return P_l_plus_1;
}
