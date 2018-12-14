from boostksp import libksp
import psutil
import os
import matplotlib.pyplot as plt
import numpy as np


source_ind = 0
sink_ind = 1

nodes = {'a': source_ind,
         'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7, 'z': sink_ind}

inv_nodes = {v: k for k, v in nodes.items()}


def make_graph(nodes):

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

    return g, edges

g, edges = make_graph(nodes)
# res = g.run()

# def print_paths(res, edges, inv_nodes):
#     for k in range(len(res)):
#         print('k={}'.format(k))
#         print('--------------')
#         for e in res[k]:
#             print('({},{})'.format(inv_nodes[edges[e][0]],
#                                 inv_nodes[edges[e][1]]))

# print_paths(res, edges, inv_nodes)

process = psutil.Process(os.getpid())
mems = []

num_iters = 50000

# for i in range(num_iters):
while True:

    # if(i == num_iters // 2):
    #     g, edges = make_graph(nodes)
    #     g.print_all()

    mem = process.memory_percent()
    # print((i, mem))
    # mems.append((i, mem))

    idx_edge_to_remove = np.random.choice(np.arange(0, len(edges)))
    edge_to_remove = edges[idx_edge_to_remove]

    g.remove_edge(*edge_to_remove[0:2])
    g.add_edge(*edge_to_remove)
    g.run()

# print_paths(res, edges, inv_nodes)
local_ram = 19
arr_mem = np.asarray([(i, t[-1]/100*local_ram*1000) for i, t in enumerate(mems)])
arr_mem[:, 1] = arr_mem[:, 1] - np.min(arr_mem[:, 1])
plt.plot(arr_mem[:, 0], arr_mem[:, 1], 'bo-')
plt.xlabel('iteration')
plt.ylabel('RAM occupancy [MB]')
plt.show()
