from boostksp import libksp


source_ind = 0
sink_ind = 1

node_start = int(2)
nodes = {'a': source_ind,
         'b': node_start,
         'c': node_start+1,
         'd': node_start+2,
         'e': node_start+3,
         'f': node_start+4,
         'g': node_start+5,
         'z': sink_ind}

inv_nodes = {v: k for k, v in nodes.items()}


g = libksp.ksp()
g.config(0,
         1,
         loglevel='info',
         l_max=10,
         min_cost=False,
         return_edges=True)

edges = list()
id_e = int(0)
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

# while True:
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

print('will remove edges...')
print(g.out_edges(v))
print('-'*50)
print('num. vertices before: {}'.format(g.num_vertices()))
print('num. edges before: {}'.format(g.num_edges()))
g.print_all()
[g.remove_edge(e[0], e[1]) for e in g.out_edges(v)]
print('num. vertices after: {}'.format(g.num_vertices()))
print('num. edges after: {}'.format(g.num_edges()))
g.print_all()
