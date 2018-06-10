from boostksp import libksp


source_ind = 0
sink_ind = 1

nodes = {'s': source_ind,
         'a': 2, 'b': 3, 'c': 4, 'd': 5, 't': sink_ind}

g = libksp.ksp.create()
g.new_graph(len(nodes.keys()))

id_e = 0
g.add_edge(nodes['s'], nodes['a'], 3, id_e)
id_e += 1
g.add_edge(nodes['s'], nodes['b'], 2, id_e)
id_e += 1
g.add_edge(nodes['a'], nodes['c'], 4, id_e)
id_e += 1
g.add_edge(nodes['b'], nodes['d'], 3, id_e)
id_e += 1
g.add_edge(nodes['c'], nodes['t'], 1, id_e)
id_e += 1
g.add_edge(nodes['d'], nodes['t'], 2, id_e)
id_e += 1
g.add_edge(nodes['b'], nodes['a'], 1, id_e)
id_e += 1
g.add_edge(nodes['c'], nodes['d'], 2, id_e)
id_e += 1
g.add_edge(nodes['b'], nodes['c'], 2, id_e)

g.set_source_id(source_ind)
g.set_sink_id(sink_ind)

g.print_edges()
g.do_ksp()
