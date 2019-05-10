#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble graph [-h|--help]

    Options:
        -h --help                       show this
"""

from docopt import docopt

from itertools import combinations, chain, zip_longest
from collections import Counter
import networkx as nx
from string import ascii_lowercase


# https://stackoverflow.com/questions/10035752/elegant-python-code-for-integer-partitioning

def coalesce(haplotype_sets):
    coalesced_haplotypes = []
    print("[+] coalescing %r" % haplotype_sets)
    for haplotype_set in haplotype_sets:
        haplotype_list = list(set(haplotype_set))
        for d_idx, a_idx in combinations(range(len(haplotype_list)), 2):
            coalesced_haplotype = [haplotype_list[a_idx] | haplotype_list[d_idx]]
            for o_idx in range(len(haplotype_list)):
                if not o_idx == a_idx and not o_idx == d_idx:
                    coalesced_haplotype.append(haplotype_list[o_idx])
            coalesced_haplotypes.append(frozenset(coalesced_haplotype))
    return set(coalesced_haplotypes)

def get_node_sample_id(l, d):
    if isinstance(l, list):
        print("[P]", l)
        return tuple([tuple([d[x] for x in _l]) for _l in l[0]])
    else:
        return frozenset([tuple([d[x] for x in _l]) for _l in l])

def get_coalescent_graph(haplotype_ids, sample_id_by_haplotype_id):
    coalescent_graph = nx.DiGraph(label0='', label1='')
    #for label in tqdm(coalesce(samples), total=formula, desc="[%] ", ncols=200):
    node_idx = 0
    coalescent_graph.add_node(node_idx, label0=haplotype_ids, label1=get_node_sample_id(haplotype_ids, sample_id_by_haplotype_id)) # tuple
    node_sample_id_counter = Counter()
    while 1:
        if len(haplotype_ids[0]) == 1:
            break
        coalesced_haplotypes = coalesce(haplotype_ids)
        for coalesced_haplotype in coalesced_haplotypes:

            node_sample_id = get_node_sample_id(coalesced_haplotype, sample_id_by_haplotype_id)
            if not node_sample_id in node_sample_id_counter:
                node_idx += 1
                coalescent_graph.add_node(node_idx, label0=coalesced_haplotype, label1=node_sample_id)
            node_sample_id_counter[node_sample_id] += 1
        haplotype_ids = list(set(coalesced_haplotypes))
    print(node_sample_id_counter)
    for node, label in coalescent_graph.nodes(data=True):
        print("[G]", node, len(label), label)

    #label_paths
    #get_all_paths

def get_samples_and_haplotypes(N, ploidy):
    #sample_ids = [letter for n, letter in zip(range(N), ascii_lowercase)] if N < len(ascii_lowercase) else [str(s) for s in range(N)]
    sample_ids = [letter for n, letter in zip(range(N), ascii_lowercase)] # if N < len(ascii_lowercase) # else ["%s.%s" % (s, n) for s in range(N) for n in range(ploidy)]
    print("[#] sample_ids", sample_ids)
    haplotype_ids = [str(h) for h in range(N * ploidy)]
    print("[#] haplotypes", haplotype_ids)
    haplotype_ids = [frozenset(str(h)) for h in range(N * ploidy)]
    sample_id_by_haplotype_id = {}
    idx = 0
    sample_ids_reversed = sample_ids[::-1]
    while 1:
        sample_id = sample_ids_reversed.pop()
        for p in range(ploidy):
           sample_id_by_haplotype_id[list(haplotype_ids[idx])[0]] = "%s.%s" % (sample_id, p + 1)
           idx += 1 
        if len(sample_ids_reversed) == 0:
            break
    print("[#] sample_id_by_haplotype_id", sample_id_by_haplotype_id)
    return sample_ids, [haplotype_ids], sample_id_by_haplotype_id

class CoalescentNode(object):
    def __init__(self, node_idx, node_haplotype):
        self.node_idx = node_idx
        self.node_haplotype = node_haplotype # list of tuples
        self.frozenset = frozenset(node_haplotype)
        self.length = 0 
        # node_haplotype
        # 0 : frozenset({(('a',), ('a',), ('b',), ('b',))} => 4
        # 1 : frozenset({('b',), ('a', 'a'), ('b',)})} => 3
        # 2 : frozenset({('a',), ('a',), ('b', 'b')})} => 3
        # 3 : frozenset({('a',), ('b',), ('a', 'b')})} => 3
        # 4 : frozenset({('a',), ('b',), ('b', 'a')})} => 3
        # 5 : frozenset({('b', 'a', 'a'), ('b',)})} => 2
        # 6 : frozenset({('a',), ('b', 'b', 'a')})} => 2
        # 7 : frozenset({('a', 'a'), ('b', 'b')})} => 2
        # 8 : frozenset({('a',), ('b', 'a', 'b')})} => 2
        # 9 : frozenset({('a', 'b'), ('b', 'a')})} => 2
        # 10: frozenset({('b', 'b'), ('a', 'a')} => 2
        # 11 : frozenset({('b', 'a', 'a', 'b')})} => 1

def make_graph():
    ploidy = 2
    N = 2
    sample_ids, haplotype_ids, sample_id_by_haplotype_id = get_samples_and_haplotypes(N, ploidy)
    coalescent_graph = get_coalescent_graph(haplotype_ids, sample_id_by_haplotype_id)
    print(coalescent_graph)
    
def main():
    args = docopt(__doc__)
    print(args)
    make_graph()

if __name__ == '__main__':
    main()