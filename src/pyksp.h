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
#include "boost_ksp.h"
#include "globals.h"

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace boost;


struct ksp {

    MyGraph * G;
    MyGraph * G_c; //with transformed edge costs
    MyGraph * G_l; //with transformed edge costs

    int n_vertices;
    Vertex source_vertex;
    Vertex sink_vertex;


    static shared_ptr<ksp> create();
    void new_graph(int n_vertices);

    bp::list do_ksp();
    Vertex add_vertex(int id, std::string str);

    void set_loglevel(unsigned int a_log_level);
    bool add_edge(int n0,
                     int n1,
                     double w,
                     int id=-1,
                     std::string str_0="",
                     std::string str_1="",
                     int label=1);
    EdgeSets augment(EdgeSets P_l,
                     EdgeSet p_inter,
                     MyGraph & g,
                     MyGraph & g_c,
                     MyGraph & g_l);
    void set_source(int id, std::string str);
    void set_sink(int id, std::string str);
    ShortestPathRes bellman_ford_shortest_paths(const MyGraph & g);
    ShortestPathRes dijkstra_shortest_paths(const MyGraph & g, int sink_id);

    void cost_transform(const std::vector<double> & distance,
                        const MyGraph & g_in, MyGraph & g_out);

    std::string hello() { return "Just nod if you can hear me!"; }

    void run(int n0,
                    int n1,
                    double w);
};

#endif
