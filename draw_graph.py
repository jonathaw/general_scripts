#!/usr/bin/env python3.5
import matplotlib.pyplot as plt
import networkx as nx
__author__ = 'jonathan'


def differ(n1, n2):
    count_1, count_2 = 0, 0
    for i, a in enumerate(n1[:3]):
        count_1 += 1 if (a == 'N' and n2[i] == 'P') or (a == 'P' and n2[i] == 'N') else 0
    for i, a in enumerate(n1[4:]):
        count_2 += 1 if (a == 'N' and n2[i] == 'P') or (a == 'P' and n2[i] == 'N') else 0
    return count_1 >= 2 and count_2 >= 1
fig = plt.gcf()
G = nx.Graph()

letters = ['N', 'P']
nodes = ['%s%s%s-%s%s' % (a, b, c, d, e) for a in letters for b in letters for c in letters
         for d in letters for e in letters]
[G.add_node(a, label=str(a)) for a in nodes]
[G.add_edge(a, b) for a in nodes for b in nodes if differ(a, b)]

all_cliques = [a for a in nx.find_cliques(G)]
max_size = max([len(a) for a in all_cliques])
max_clique = [clq for clq in all_cliques if len(clq) == max_size][0]

clique = [(u, v) for (u, v) in G.edges() if u in max_clique and v in max_clique]

pos=nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, c='b', alpha=0.4)
nx.draw_networkx_nodes(G, pos, nodelist=max_clique, node_color='cornflowerblue', alpha=0.8)

nx.draw_networkx_edges(G, pos, edge_color='gray', alpha=0.2)
nx.draw_networkx_edges(G, pos, edgelist=clique, edge_color='black', alpha=0.9, width=2.0)

nx.draw_networkx_labels(G, pos, c='r', font_size=20, font_family='sans-serif')
fig.set_size_inches(18.5, 10.5)
plt.savefig('/home/labs/fleishman/jonathaw/graph.png', dpi=100)

plt.show()
#!/usr/bin/env python
def main():
    import matplotlib.pyplot as plt
    import argparse
    import numpy as np
    parser = argparse.ArgumentParser()
    parser.add_argument('-eq')
    parser.add_argument('-rngx', nargs=2, default=[-10, 10])
    parser.add_argument('-rngy', nargs=2, default=[-10, 10])
    args = vars(parser.parse_args())
    for i, t in enumerate(args['rngx']):
        if str(t)[0] == 'm':
            args['rngx'][i] = - int(t[1:])
        else:
            args['rngx'][i] = int(t)
    for i, t in enumerate(args['rngy']):
        if str(t)[0] == 'm':
            args['rngy'][i] = - int(t[1:])
        else:
            args['rngy'][i] = int(t)
    print args
    x = np.array(np.arange(args['rngx'][0], args['rngx'][1], 0.01))
    print x
    res = eval(args['eq'])
    print 'res', res
    plt.xlim(args['rngx'])
    plt.ylim(args['rngy'])
    plt.plot(x, res)
    plt.show()


if __name__ == '__main__':
    main()
