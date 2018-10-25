from boostksp import libksp


source_ind = 0
sink_ind = 1

nodes = {'a': source_ind,
         'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7, 'z': sink_ind}

inv_nodes = {v: k for k, v in nodes.items()}


g = libksp.ksp()
g.config(0,
         1,
         loglevel='trace',
         min_cost=False,
         return_edges=True)

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

g.set_source(source_ind, 'y')
g.set_sink(sink_ind, 'z')

for e in edges:
    g.add_edge(*e)

res = g.run()


def print_paths(res, edges, inv_nodes):
    for k in range(len(res)):
        print('k={}'.format(k))
        print('--------------')
        for e in res[k]:
            print('({},{})'.format(inv_nodes[edges[e][0]],
                                inv_nodes[edges[e][1]]))

print_paths(res, edges, inv_nodes)

v = source_ind
print('out_edges from vertex {}'.format(v))
print(g.out_edges(v))

print('testing removal of vertex...')
print('-'*50)
print('num. vertices before: {}'.format(g.num_vertices()))
print('num. edges before: {}'.format(g.num_edges()))
g.remove_vertex(source_ind)
print('num. vertices after: {}'.format(g.num_vertices()))
print('num. edges after: {}'.format(g.num_edges()))
