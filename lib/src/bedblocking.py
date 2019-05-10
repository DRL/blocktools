from pandas import read_csv
from numpy import int as npint

from os.path import join, dirname, realpath, isdir, isfile
from os import getcwd
from collections import defaultdict, OrderedDict
from sys import exit
from shutil import rmtree
import tqdm as tqdm
from lib.functions import create
#from graph_tool.all import *
import networkx as nx

def evaluate_edges():
    # has to return weight of edge
    #   - pop_weight = set of pops
    # construct a dict of frozensets of all possible traversals based on decrease in sample_sets until min_samples/min_pops

print("Nodes in G: ", G.nodes(data=True))
print("Edges in G: ", G.edges(data=True))


from itertools import combinations, product
from functools import reduce  
import operator
from pprint import pprint

def getFromDict(dataDict, sample_key):
    return reduce(operator.getitem, sample_key, dataDict)

def setInDict(dataDict, sample_key, value):
    getFromDict(dataDict, sample_key[:-1])[sample_key[-1]] = value

sample_ids = set(['A', 'B', 'C', 'X', 'Y', 'Z'])
min_sample = 3
d = {}
sample_keys = []
for sample_count in range(len(sample_ids), min_sample-1, -1):
    sample_sets = [frozenset(x) for x in combinations(sample_ids, sample_count)]
    print("#", sample_sets)
    for sample_set in sample_sets:
        sample_keys.append([sample_set])
        if sample_set < sample_key[-1]:
            sample_key.append(sample_set)
        if sample_set < sample_key[0]:
            sample_keys.append([sample_key[0], sample_set])
        print("M", sample_keys)
        
            
        
            

class BedObj(object):
    #__slots__ = ['chrom', 'start', 'end', 'pair_idxs', 'length']
    
    def __init__(self, chrom, start, end, sample_idxs, pair_idxs, length):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.length = int(length)
        self.pair_idxs = set(pair_idxs)
        self.sample_idxs = set(sample_idxs)

class BedGraph(object):
    __slots__ = ["contig_id", "bedObjs", "idx"]
    def __init__(self, seq_id, bedObjs):
        self.seq_id = seq_id
        self.bedObjs = bedObjs # list

        self.graph = nx.DiGraph()

    def build_graph(self):
        
        for idx, bedObj in enumerate(bedObjs):
            v = self.graph.add_vertex(idx)
            vprop_length[v] = bedObj.length
            vprop_valid[v] = True

    def get_best_path(self, sample_counter, pop_counter):
        if all([sample_counter, pop_counter]):
            pass
        elif sample_counter:
            pass
        elif pop_counter:
            pass
        else:
            # first implement without frequency filters 

            pass

class ParameterObj(object):
    def __init__(self, args):
        self.args = args
        self.data_path = join(dirname(realpath(__file__)), "../data/")
        self.plot_size = (24, 12)
        self.plot_font_size = 18
        self.out_dir = join(getcwd(), args['--out_dir'])
        self.config_dir = join(self.out_dir, "0_config")

        self.config_f = self.check_parameter(name='--config_file', param_type='file')
        self.multibed_f = self.check_parameter(name='--config_file', param_type='file')
        self.min_samples = self.check_parameter(name="--min_samples", param_type='posint')
        self.min_pops = self.check_parameter(name="--min_pops", param_type='posint')

        self.sample_ids_by_pop_id = defaultdict(set)
        self.sample_ids = set()
        self.pop_ids_by_sample_id = {}
        self.pop_ids = set()
        self.length_by_sequence_id = OrderedDict()
        
    def parse_samples_f(self):
        sample_df = read_csv( \
            self.samples_f, \
            sep=",", \
            usecols=[0, 1], \
            names=['sample_id', 'pop_id'], \
            dtype={'sample_id': 'category', 'pop_id': 'category'})
        for sample_id, pop_id in tqdm(sample_df.values.tolist(), total=len(sample_df.index), desc="[%] ", ncols=200):
            self.sample_ids_by_pop_id[pop_id].add(sample_id)
            self.pop_id_by_sample_id[sample_id] = pop_id
            self.sample_ids.add(sample_id)
            self.pop_ids.add(pop_id)
        if len(self.sample_ids) < self.min_samples:
            exit("[X] %r contains fewer samples than '--min_samples' (%r) " % (self.samples_f, self.min_samples))
        if len(self.pop_ids) < self.min_pops:
            exit("[X] %r contains fewer populations than '--min_pops' (%r) " % (self.samples_f, self.min_pops))

    def parse_multibed_f(self):
        df = read_csv( \
            self.multibed_f, \
            sep="\t", \
            usecols=[0, 1, 2, 4], \
            names=['chrom', 'start', 'end', 'samples'], \
            skiprows=1, \
            header=None, \
            dtype={ \
                'chrom': 'category', \
                'start': npint, \
                'end': npint, \
                'samples': 'category'\
            })
        df = df[df['chrom'].isin(self.length_by_sequence_id)]
        df['length'] =  df['end'] - df['start']
        # df = df[df['length'] >= parameterObj.min_interval_len]
        df['samples_idxs'] = df['samples'].apply(\
            generate_sample_idxs, \
            sample_idx_by_sample_id=parameterObj.sample_idx_by_sample_id\
            )
        df['pair_idxs'] = df['samples'].apply(\
            generate_pair_idxs, \
            pops_count=parameterObj.pops_count, \
            pair_idx_by_pair_ids=parameterObj.pair_idx_by_pair_ids\
            )
        # Drop intervals that don't affect pairs
        # df = df.dropna()
        df['distance'] = numpy.where((df['chrom'] == df['chrom'].shift(-1)), df['start'].shift(-1) - df['end'], numpy.nan)
        for chrom, start, end, samples, length, sample_idxs, pair_idxs, distance in tqdm(df.values.tolist(), total=len(df.index), desc="[%] ", ncols=200):
            if not pair_idxs is numpy.nan:
                pair_count = len(pair_idxs)
                coverageObj.add_pair_region(pair_count, length)
                if length >= parameterObj.block_length:
                    coverageObj.add_pair_region_block_size(pair_count, int(length / parameterObj.block_length) * parameterObj.block_length)
                if length >= parameterObj.min_interval_len:
                    bedObj = BedObj(chrom, start, end, pair_idxs, length) 
                    if not regionBatchObj.contig_id:
                        regionBatchObj.add_bedObj_to_batch(bedObj)
                    else:
                        if chrom == regionBatchObj.contig_id and not numpy.isnan(distance) and int(distance) <= parameterObj.max_interval_distance:
                            regionBatchObj.add_bedObj_to_batch(bedObj)
                        else:
                            regionBatchObjs.append(regionBatchObj)
                            idx += 1
                            regionBatchObj = RegionBatchObj(idx)
                            regionBatchObj.add_bedObj_to_batch(bedObj)
                    #if numpy.isnan(distance) or int(distance) > parameterObj.max_interval_distance:
            
            sample_count = len(sample_idxs)
            coverageObj.add_sample_region(sample_count, length)
            if length >= parameterObj.block_length:
                coverageObj.add_sample_region_block_size(sample_count, int(length / parameterObj.block_length) * parameterObj.block_length)
        if regionBatchObj.contig_id:
            regionBatchObjs.append(regionBatchObj)
        #for x in regionBatchObjs:
        #    print(x)
            #for y in x.bedObjs:
            #    print(y)
    
        # awk '{if($3-$2>=65){a[$4]+=$3-$2;}}END{for (i in a)print i", "a[i];}' data/input/tiny/hmel.autosomes.tiny.multiinter.samples_as_string.bed
        return regionBatchObjs, coverageObj
                
    def create_out_dir(self):
        if isdir(self.out_dir):
            print("[!] Directory %r exists. Deleting directory ..." % self.out_dir)
            rmtree(self.out_dir)
        create(out_f=None, path=self.out_dir, header=None, lines=None)

    def check_out_dir(self):
        if isdir(self.out_dir):
            return True
        return False

    def check_parameter(self, name='', param_type=''):
        value = self.args.get(name, None)
        if value:
            if param_type == 'bool':
                if value is True or value is False:
                    return value
            if param_type == "posint":
                try:
                    param = int(value)
                    if param > 0:
                        return param
                except TypeError:
                    pass
                exit("[X] %r must be a positive integer. Was %r." % (name, value))
            elif param_type == "int":
                try:
                    return int(value)
                except TypeError:
                    pass
                exit("[X] %r must be an integer. Was %r." % (name, value))
            elif param_type == "file":
                return self.check_file(value)
            elif param_type == "dir":
                try: 
                    if isdir(value):
                        return value
                except TypeError:
                    pass
                exit("[X] Directory not found: %r." % value)
            elif param_type == "fraction":
                try:
                    param = float(value)
                    if 0 <= param <= 1:
                        return param
                except TypeError:
                    pass
                exit("[X] %r must be a fraction. Was %r." % (name, value))
            elif param_type == 'plotfmt':
                plot_fmts = ['png', 'pdf', 'svg']
                if value in set(plot_fmts):
                    return value
                exit("[X] %r must be one of: %r. Was %r." % (name, ", ".join(plot_fmts), value))
        return None

    def check_file(self, infile):
        try: 
            if isfile(infile):
                return infile
        except TypeError:
            pass
        exit("[X] File not found: %r." % infile)