from boostksp import libksp


source_ind = 0
sink_ind = 1

nodes = {'s': source_ind,
         'b': 2, 'c': 3, 'd': 4, 'e': 5, 't': sink_ind}

g = libksp.ksp.create()
g.new_graph(len(nodes.keys()))

g.add_edge(nodes['s'], nodes['b'], 4)
g.add_edge(nodes['s'], nodes['c'], 2)
g.add_edge(nodes['b'], nodes['c'], 5)
g.add_edge(nodes['c'], nodes['e'], 3)
g.add_edge(nodes['b'], nodes['d'], 10)
g.add_edge(nodes['e'], nodes['d'], 4)
g.add_edge(nodes['d'], nodes['t'], 11)

g.set_source_id(source_ind)
g.set_sink_id(sink_ind)
g.bellman_ford_shortest_paths()
