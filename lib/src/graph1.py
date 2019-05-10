#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble graph -n INT -p INT [-h|--help]

    Options:
        -h --help                       show this
        -n, --samples INT               Number of samples
        -p, --ploidy INT                Ploidy of samples
"""

from docopt import docopt

from itertools import combinations
from collections import Counter
from networkx.drawing.nx_agraph import graphviz_layout
import networkx as nx
from string import ascii_lowercase
import matplotlib as mat
mat.use("agg")
import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/10035752/elegant-python-code-for-integer-partitioning

def coalesce_node(parent_node):
    child_nodes = []
    for d_idx, a_idx in combinations(range(len(parent_node.haplotypes)), 2):
        child_haplotypes = [haplotype for idx, haplotype in enumerate(parent_node.haplotypes) if not idx == d_idx and not idx == a_idx]
        child_haplotypes.append("".join(sorted(parent_node.haplotypes[a_idx] + parent_node.haplotypes[d_idx])))
        child_node = Node(idx=None, haplotypes=child_haplotypes, parent=parent_node.idx)
        child_nodes.append(child_node)
    return Counter(child_nodes)

class Node(object):
    def __init__(self, haplotypes='', idx=None, parent=None):
        self.idx = idx
        self.haplotypes = tuple(sorted(haplotypes))
        self.parent = parent

    def __eq__(self, other):
        if isinstance(other, Node):
            return self.haplotypes == other.haplotypes

    def __str__(self):
        return "idx=%s haplotypes=%s parent=%s" % (self.idx, self.haplotypes, self.parent)

    def __hash__(self):
        return hash(self.haplotypes)

def build_state_graph(source_node):
    
    state_graph = nx.DiGraph(label='')
    source_node.idx = 0
    state_graph.add_node(source_node.idx, label=source_node.haplotypes)
    node_queue = [source_node]
    nodes_created = set() 
    edges = []
    node_idx = 0
    node_idx_by_haplotype = {}
    while 1:
        if len(node_queue) == 0:
            break
        parent_node = node_queue.pop(0)
        child_node_counter = coalesce_node(parent_node)
        node_queue += list(child_node_counter.keys())
        for child_node, count in child_node_counter.items():
            child_node.idx = node_idx_by_haplotype.get(child_node.haplotypes, None)
            if child_node.idx == None:
                node_idx += 1
                child_node.idx = node_idx
                state_graph.add_node(child_node.idx, label=child_node.haplotypes)
                node_idx_by_haplotype[child_node.haplotypes] = child_node.idx
            edges.append((child_node.parent, node_idx_by_haplotype[child_node.haplotypes], {'count': count}))
            nodes_created.add(child_node)
    state_graph.add_edges_from(edges)
    return state_graph 

def get_source_node(N, ploidy):
    sample_ids = [letter for n, letter in zip(range(N), ascii_lowercase)] 
    print("[#] sample_ids", sample_ids)
    haplotypes = sorted([s_id for s_id in sample_ids for p in range(ploidy)])
    print("[#] haplotypes", haplotypes)
    source_node = Node(idx=0, haplotypes=haplotypes)
    #print('[0]', source_node)
    return source_node

def plot_state_graph(state_graph):
    fig = plt.figure(figsize=(14,14), dpi=400, frameon=False)
    options = {
        'node_color': 'C0',
        'node_size': 700,
        'with_labels': True,
    }
    pos=graphviz_layout(state_graph, prog='dot')
    node_labels = {n: "[%s] %s" % (n, d.get('label', 'None')) for n, d in state_graph.nodes(data=True)}
    nx.draw_networkx_nodes(state_graph, pos, **options)
    nx.draw_networkx_labels(state_graph, pos, labels={n: lab for n, lab in node_labels.items() if n in pos}, font_size=16, font_family='sans-serif')
    edge_labels = {(source, sink): d.get('count', 'None') for source, sink, d in state_graph.edges(data=True)}
    nx.draw_networkx_edges(state_graph, pos, width=2)
    nx.draw_networkx_edge_labels(state_graph, pos, edge_labels=edge_labels, label_pos=0.5, font_size=20)
    out_f = "test.png"
    plt.axis('off')
    fig.savefig(out_f, format="png")
    plt.close(fig)


def make_graph(args):
    ploidy = int(args['--ploidy'])
    N = int(args['--samples'])
    print("[ARGS] N=%s, p=%s" % (N, ploidy))
    source_node = get_source_node(N, ploidy)
    state_graph = build_state_graph(source_node)
    print("# [GRAPH]")
    for idx, node in state_graph.nodes(data=True):
        print("[%s] %s" % (idx, node['label']))
    plot_state_graph(state_graph)
    print("# [PATHS]")
    for idx, path in enumerate(nx.all_simple_paths(state_graph, source=0, target=len(state_graph)-1)):
        print("[%s] %s weight=%s" % (idx, path, sum([state_graph[x][y]['count'] for x, y in pairwise(path)])))

def pairwise(iterable):
    it = iter(iterable)
    a = next(it, None)
    for b in it:
        yield (a, b)
        a = b

def main():
    args = docopt(__doc__)
    make_graph(args)

if __name__ == '__main__':
    main()