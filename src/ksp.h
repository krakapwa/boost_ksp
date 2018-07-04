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

#ifndef PYKSP_H
#define PYKSP_H

#include <iostream>
#include <limits>
#include "globals.h"
#include "utils.h"

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace boost;


class Ksp {

    public:
        Ksp();
        void config(int source_vertex_id,
                    int sink_vertex_id,
                    int l_max,
                    std::string source_vertex_name,
                    std::string sink_vertex_name,
                    std::string loglevel,
                    bool min_cost);
        Vertex add_vertex(int id, std::string str);

        void set_loglevel(std::string);
        bool add_edge(int n0,
                            int n1,
                            double w,
                            int id=-1,
                            std::string str_0="",
                            std::string str_1="",
                            int label=1);
        void remove_edge(int u, int v);

        void set_source(int id, std::string str);
        void set_sink(int id, std::string str);
        void set_label_all_edges(int label);
        int num_vertices();
        int num_edges();


        bp::list run();

    private:

        MyGraph * G;
        MyGraph * G_c; //with transformed edge costs

        int n_vertices;
        Vertex source_vertex;
        Vertex sink_vertex;

        int l_max;

        bool min_cost;
        bool return_edges;

        double cost;
        double new_cost; // will store two consecutives costs for comparison

        void new_graph();


        std::tuple<EdgeSet, bool, std::vector<double>>
        bellman_ford_shortest_paths(const MyGraph & g);

        std::tuple<EdgeSet, bool, std::vector<double>>
        dijkstra_shortest_paths(const MyGraph & g, Vertex source_vertex);

        void cost_transform(const std::vector<double> & distance,
                            MyGraph & g_out);

        std::string hello() { return "Just nod if you can hear me!"; }

};

#endif
