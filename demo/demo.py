from boostksp import libksp


source_ind = 0
sink_ind = 1

nodes = {'a': source_ind,
         'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7, 'z': sink_ind}

inv_nodes = {v: k for k, v in nodes.items()}

g = libksp.ksp()
g.config(source_vertex_id=0,
         sink_vertex_id=1,
         loglevel='debug')

edges = list()
id_e = 0
edges.append((nodes['a'], nodes['b'], 1, id_e, 'a', 'b'))
id_e += 1
edges.append((nodes['b'], nodes['c'], 1, id_e, 'b', 'c'))
id_e += 1
edges.append((nodes['c'], nodes['d'], 1, id_e, 'c', 'd'))
id_e += 1
edges.append((nodes['d'], nodes['z'], 1, id_e, 'd', 'z'))
id_e += 1
edges.append((nodes['a'], nodes['e'], 1, id_e, 'a', 'e'))
id_e += 1
edges.append((nodes['e'], nodes['f'], 3, id_e, 'e', 'f'))
id_e += 1
edges.append((nodes['f'], nodes['z'], 5, id_e, 'f', 'z'))
id_e += 1
edges.append((nodes['e'], nodes['b'], 1, id_e, 'e', 'b'))
id_e += 1
edges.append((nodes['b'], nodes['f'], 1, id_e, 'b', 'f'))
id_e += 1
edges.append((nodes['f'], nodes['d'], 2, id_e, 'f', 'd'))
id_e += 1
edges.append((nodes['a'], nodes['g'], 7, id_e, 'a', 'g'))
id_e += 1
edges.append((nodes['g'], nodes['z'], 2, id_e, 'g', 'z'))
id_e += 1
edges.append((nodes['c'], nodes['g'], 1, id_e, 'c', 'g'))
id_e += 1

g.set_source(source_ind, 'a')
g.set_sink(sink_ind, 'z')

for e in edges:
    g.add_edge(*e)

res = g.run()

for k in range(len(res)):
    print('k={}'.format(k))
    print('--------------')
    for e in res[k]:
        print('({},{})'.format(inv_nodes[edges[e][0]],
                               inv_nodes[edges[e][1]]))
