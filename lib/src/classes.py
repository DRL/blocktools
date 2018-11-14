from itertools import product
from collections import OrderedDict, defaultdict, deque, Counter
import operator
from tqdm import tqdm
from lib.src.code import *

class ProfileObj():
    '''
    Could be replaced with a NamedTuple
    '''
    __slots__ = ['fixed', 'hetA', 'hetAB', 'hetB', 'invariant', 'missing', 'multiallelic']

    def __init__(self, profile_tuple):
        self.fixed = profile_tuple[0]
        self.hetA = profile_tuple[1]
        self.hetAB = profile_tuple[2]
        self.hetB = profile_tuple[3]
        self.invariant = profile_tuple[4]
        self.missing = profile_tuple[5]
        self.multiallelic = profile_tuple[6]

    def mutuple(self):
        return (self.fixed, self.hetA, self.hetAB, self.hetB)

    def tuple(self):
        return (self.fixed, self.hetA, self.hetAB, self.hetB, self.invariant, self.missing, self.multiallelic)

    def __str__(self):
        return str(self.tuple())

    def __add__(self, otherProfileObj):
        '''
        python -m timeit "x=(1,2,3,4,5); y=(0,1,2,3,4); tuple(map(sum, zip(x, y)))"
            1000000 loops, best of 3: 0.923 usec per loop
        python -m timeit "x=(1,2,3,4,5); y=(0,1,2,3,4); tuple(_x + _y for _x, _y in zip(x, y))"
            1000000 loops, best of 3: 0.943 usec per loop
        python -m timeit "import operator; x=(1,2,3,4,5); y=(0,1,2,3,4); tuple(map(operator.add, x, y))"
            1000000 loops, best of 3: 0.783 usec per loop 
        '''
        return ProfileObj(tuple(map(operator.add, self.tuple(), otherProfileObj.tuple())))

    def __radd__(self, otherProfileObj):
        if isinstance(otherProfileObj, int):
            return self
        else:
            return self.__add__(otherProfileObj)

    def __truediv__(self, integer):
        return ProfileObj(tuple([x/integer for x in self.tuple()]))


class CoordinateTransformObj(object):
    def __init__(self):
        self.tuples_by_seq_id = defaultdict(list)
    
    def add_tuple(self, seq_id, seq_start, seq_end, seq_orientation, chrom, chrom_start, chrom_end):
        self.tuples_by_seq_id[seq_id].append((seq_id, seq_start, seq_end, seq_orientation, chrom, chrom_start, chrom_end))

    def transform_coordinates(self, seq_id, seq_start, seq_end):
        if seq_id in self.tuples_by_seq_id:
            for _tuple in self.tuples_by_seq_id[seq_id]:
                if seq_start >= _tuple[1] and seq_end <= _tuple[2]:
                    if _tuple[3] == "+":
                        return _tuple[4], _tuple[5] + seq_start, _tuple[5] + seq_end
                    else:
                        '''
                        GFF :   chr15   5188091 5243216 Hmel215037  1   55126   -
                        BED :   chr15   5188090 5243216 Hmel215037  0   55126   -
                        Simon:  Hmel215037 549 672 + |  ('chr15', 5242545, 5242668) 
                        Dom:    Hmel215037 549 672   |  ('chr15', 5242544, 5242667)
                        '''
                        return _tuple[4], _tuple[6] - seq_end, _tuple[6] - seq_start
        return None, None, None


class ParameterObj(object):
    '''
    Class containing all parameters necessary for:
        - blocks
        - variants
        - windows
        - plots
        - picking up after each checkpoint (!!!)
    '''
    def __init__(self, args):
        # input files
        self.genome_f = check_file(args['--genome'])
        self.multibed_f = check_file(args['--multibed'])
        self.vcf_f = check_file(args['--vcf'])
        # analysis parameters
        self.threads = int(args['--threads'])
        self.algorithm = args['--algorithm']
        self.block_length = int(args['--block_length'])
        self.min_interval_len = int(args['--min_interval_len'])
        self.max_interval_distance = int(args['--max_interval_distance'])
        self.block_span = self.max_interval_distance

        #self.modifier = float(args['--modifier'])
        self.window_size = int(args['--window_size']) if '--window_size' in args else None
        self.window_overlap = int(args['--overlap']) if '--overlap' in args else None
        # Samples/Pops
        self.sample_ids_by_population = args['--populations']
        self.pops_count = len(self.sample_ids_by_population)
        if not self.pops_count == 2:
            raise PopsCountException()
        self.sample_ids = [sample_id for sample_ids in sorted(self.sample_ids_by_population.values()) for sample_id in sample_ids]
        self.samples_count = len(self.sample_ids)
        self.sample_id_by_sample_idx = {idx: sample_id for idx, sample_id in enumerate(self.sample_ids)}
        self.sample_idx_by_sample_id = {sample_id: idx for idx, sample_id in enumerate(self.sample_ids)}
        self.pair_ids = [(x) for x in product(*self.sample_ids_by_population.values())]
        self.pairs_count = len(self.pair_ids)
        #self.intersection_threshold = self.modifier * self.pairs_count
        self.pair_idxs = [pair_idx for pair_idx, pair_id in enumerate(self.pair_ids)]
        self.pair_idx_by_pair_ids = {frozenset(pair_id): pair_idx for pair_idx, pair_id in zip(self.pair_idxs, self.pair_ids)}
        self.pair_ids_by_pair_idx = {pair_idx: pair_id for pair_idx, pair_id in zip(self.pair_idxs, self.pair_ids)}
        self.sample_idxs_by_pair_idx = {pair_idx: (self.sample_idx_by_sample_id[pair_id[0]], self.sample_idx_by_sample_id[pair_id[1]]) for pair_idx, pair_id in zip(self.pair_idxs, self.pair_ids)}
        # Output files
        self.outprefix = args['--outprefix']
        self.block_matrix_yaml_f = "%s.block_matrix.yaml" % (self.outprefix)
        self.block_bed_f = "%s.block.bed" % (self.outprefix)
        self.block_tsv_f = "%s.block.tsv" % (self.outprefix)
        self.profiles_by_block_id_tsv_f = "%s.profiles_by_block.tsv" % (self.outprefix)
        self.profiles_summary_tsv_f = "%s.profiles_summary.tsv" % (self.outprefix)
        self.window_metrics_tsv_f = "%s.window_metrics.tsv" % (self.outprefix)
        self.window_bed_f = "%s.window.bed" % (self.outprefix)
        self.window_profiles_tsv_f = "%s.window.profiles.tsv" % (self.outprefix)

class RegionBatchObj(object):
    __slots__ = ["contig_id", "bedObjs", "idx"]

    def __init__(self, idx):
        self.contig_id = None
        self.bedObjs = deque()
        self.idx = idx

    def __str__(self):
        return "[I] : c=%s\tbedObjs=%s\tlength=%s\tspan=%s" % (self.contig_id, len(self.bedObjs), self.length(), self.span())

    def add_bedObj_to_batch(self, bedObj):
        if self.contig_id == None:
            self.contig_id = bedObj.chrom
        self.bedObjs.append(bedObj)

    def length(self):
        try:
            return sum([len(bedObj) for bedObj in self.bedObjs])
        except TypeError:
            return 0

    def span(self):
        try:
            return (self.bedObjs[-1].end - self.bedObjs[0].start)
        except IndexError:
            return 0

class BedObj(object):
    __slots__ = ['chrom', 'start', 'end', 'pair_idxs', 'length']
    
    def __init__(self, chrom, start, end, pair_idxs, length):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.length = int(length)
        self.pair_idxs = set(pair_idxs)

    def __str__(self):
        return "\t".join([self.chrom, str(self.start), str(self.end), str(self.length), str(self.pair_idxs)]) 

class BlockObj(object):

    __slots__ = [\
        "contig_id", \
        "block_id", \
        "pair_idxs", \
        "sample_idxs", \
        "start", \
        "end", \
        "length", \
        "span", \
        "score", \
        "needed", \
        "bed_tuples", \
        "genotypes_by_sample_idx", \
        "profileObj_by_pair_idx", \
        ]

    def __init__(self, block_id, block_length):
        self.contig_id = block_id.split(".")[0]
        self.block_id = block_id
        self.pair_idxs = None
        self.sample_idxs = None
        self.start = None
        self.end = None
        self.length = 0
        self.span = 0 
        self.score = 0.0
        self.needed = block_length
        self.bed_tuples = [] # list of tuples (contig_id, start, end) of consecutive regions

        self.genotypes_by_sample_idx = {}  # dict of lists
        self.profileObj_by_pair_idx = {}  # dict of Counters 

    def __str__(self):
        return "[B] ID=%s %s %s LEN=%s SPAN=%s SCORE=%s [P]=%s" % (self.block_id, self.start, self.end, self.length, self.span, self.score, self.pair_idxs)

    def __nonzero__(self):
        if self.length:
            return True
        return False

    def add_bedObj(self, bedObj, parameterObj):
        '''
        Function for adding a bedObj to the blockObj

        [parameters]
            - bedObj to be added
            - parameterObj

        [returns]
            a) None (if bedObj has been consumed)
            b) bedObj
                b1) original bedObj (if span-violation)
                b2) remainder bedObj (if bedObj.length > blockObj.needed)
        
        [comments]
            - span-violation:
                if blockObj.span > parameterObj.block_span:
                    - original bedObj is returned, blockObj.score is set to 0.0
            - blockObj.needed: allows distinction between 
                a) finished block: blockObj.needed = 0
                b) virgin block: blockObj.needed = parameterObj.block_length
                c) started block: 0 < blockObj.needed < parameterObj.block_length
            - blockObj.score: 
                a) if blockObj.needed == parameterObj.block_length (virgin block):
                    blockObj.score = (bedObj.pair_idxs / parameterObj.pairs_count) * (min(bedObj.length, required_length) / parameterObj.block_length)
                b) if blockObj.needed < parameterObj.block_length (non-virgin block):
                    blockObj.score = (len(blockObj.pair_idxs.intersection(bedObj.pair_idxs)) / parameterObj.pairs_count) * (min(bedObj.length, required_length) / parameterObj.block_length)
                c) if span-violation:
                    blockObj.score = 0.0
        '''
        interval_length = min(self.needed, bedObj.length)
        _end = bedObj.start + interval_length
        try:
            _span = _end - self.start # TypeError: int() argument must be a string, a bytes-like object or a number, not 'NoneType' 
            #print("self.span:", self.span, ", _end:", _end, ", self.start:", self.start) 
            #print("_span = _end - self.start =", _span) 
        except TypeError:
            self.start = bedObj.start
            _span = _end - bedObj.start

            #print("self.span:", self.span, ", _end:", _end, ", bedObj.start:", bedObj.start) 
            #print("_span = _end - bedObj.start =", _span) 
        if _span > parameterObj.block_span:
            self.score = 0.0
            return bedObj
        # No span violation
        try:
            self.pair_idxs = self.pair_idxs.intersection(bedObj.pair_idxs) # AttributeError: 'NoneType' object has no attribute 'intersection'
        except AttributeError:
            # populate virgin block
            self.pair_idxs = bedObj.pair_idxs
            self.start = bedObj.start
        self.end = _end 
        self.span = _span
        self.length += interval_length
        self.score = (len(self.pair_idxs) / parameterObj.pairs_count) * (self.length / parameterObj.block_length)
        self.sample_idxs = pairs_to_samples(self.pair_idxs, parameterObj)
        self.needed -= interval_length
        # check whether gap: False if 0 / True if int / TypeError if block.end = None 
        try:
            if bool(self.end - bedObj.start): # TypeError: unsupported operand type(s) for -: 'NoneType' and 'int'
                # new bed_tuple
                self.bed_tuples.append((self.contig_id, bedObj.start, bedObj.end)) 
            else:
                # update last bed_tuple
                self.bed_tuples[-1] = (self.contig_id, self.start, bedObj.end) 
        except TypeError:
            # new bed_tuple
            self.bed_tuples.append((self.contig_id, bedObj.start, bedObj.end))
        if interval_length < bedObj.length:
            return BedObj(bedObj.chrom, (bedObj.start + interval_length), bedObj.end, bedObj.pair_idxs, (bedObj.length - interval_length))    
        else:
            return None
            
            

class BlockDataObj(object):
    __slots__ = ["blockObjs", "blockObj_idxs", "idx_list_by_pair_idx"]

    def __init__(self, parameterObj):
        self.blockObjs = []
        self.blockObj_idxs = [] # list of start_idxs of contigs, hacky ...
        self.idx_list_by_pair_idx = OrderedDict({pair_idx : [] for pair_idx in parameterObj.pair_idxs})

    def __len__(self):
        return len(self.blockObjs)

    def add_blockObjs(self, blockObjs):
        '''
        takes list of blockObjs
        ''' 
        blockObjs.sort(key=lambda i: (i.contig_id, i.start)) # so that they are ordered by contig/start ... 
        last_contig_id, start, end = blockObjs[0].contig_id, 0, 0
        for idx, blockObj in enumerate(blockObjs):
            if not blockObj.contig_id == last_contig_id:
                end = idx
                self.blockObj_idxs.append((start, end))
                start = idx
                last_contig_id = blockObj.contig_id
            self.blockObjs.append(blockObj)
            for pair_idx in blockObj.pair_idxs:
                self.idx_list_by_pair_idx[pair_idx].append(len(self.blockObjs) - 1) 
        
    def write_block_span_tsv(self, parameterObj):
        lines_blockmatrix_tsv = []
        lines_blockmatrix_tsv.append("#%s" % ("\t".join(["pair_idx", "pair_id", "block_count", "bases"])))
        for pair_idx, idx_list in self.idx_list_by_pair_idx.items():
            pair_id = "_vs_".join(sorted([parameterObj.pair_ids_by_pair_idx[pair_idx][0], parameterObj.pair_ids_by_pair_idx[pair_idx][1]]))
            blocks = len(idx_list)
            bases = (len(idx_list) * parameterObj.block_length)
            lines_blockmatrix_tsv.append("\t".join([str(x) for x in [pair_idx, pair_id, blocks, bases]]))
        fn_block_tsv = parameterObj.block_tsv_f
        with open(fn_block_tsv, 'w') as fh_blockmatrix_tsv:
            fh_blockmatrix_tsv.write("\n".join(lines_blockmatrix_tsv))
        return fn_block_tsv

    def write_bed(self, parameterObj):
        lines_bed = []
        lines_bed.append("#%s" % ("\t".join(["CHROM", "START", "END", "block_id", "length", "span", "count_samples", " count_pairs", "samples", "pairs"])))
        for blockObj in sorted(self.blockObjs, key=lambda i: (i.contig_id, i.start)):
            lines_bed.append("\t".join([str(x) for x in [blockObj.contig_id, blockObj.start, blockObj.end, blockObj.block_id, blockObj.length, blockObj.span, len(blockObj.sample_idxs), len(blockObj.pair_idxs), ",".join([str(x) for x in blockObj.sample_idxs]), ",".join([str(x) for x in blockObj.pair_idxs])]]))
        fn_bed = parameterObj.block_bed_f
        with open(fn_bed, 'w') as fh_bed:
            fh_bed.write("\n".join(lines_bed))
        return fn_bed

    def write_variant_output(self, parameterObj):
        '''
        BlockID -> 
        '''
        data_profiles_by_block_id = []
        data_profiles_by_block_id.append("#%s" % ("\t".join(["block_id", "pair_idx", "fixed", "hetA", "hetAB", "hetB", "invariant", "multiallelic", "missing"])))

        profileObjs_by_pair_idx = defaultdict(list)
        for blockObj in self.blockObjs:
            for pair_idx, profileObj in blockObj.profileObj_by_pair_idx.items():
                #print("\t".join([str(x) for x in profileObj.tuple()]) )
                profileObjs_by_pair_idx[pair_idx].append(profileObj)
                data_profiles_by_block_id.append("\t".join([ \
                    blockObj.block_id, \
                    str(pair_idx), \
                    "\t".join([str(x) for x in profileObj.tuple()]) \
                ]))
        
        fn_profiles_by_block_id_tsv = parameterObj.profiles_by_block_id_tsv_f
        with open(fn_profiles_by_block_id_tsv, 'w') as fh_profiles_by_block_id_tsv:
            fh_profiles_by_block_id_tsv.write("\n".join(data_profiles_by_block_id) + "\n")

        data_profiles_summary = []
        data_profiles_summary.append("#%s" % ("\t".join(["pair_idx", "blocks", "bases", "fixed", "hetA", "hetAB", "hetB", "invariant", "multiallelic", "missing"])))

        for pair_idx in parameterObj.pair_idxs:
            profileObjs = profileObjs_by_pair_idx[pair_idx]
            bases = len(profileObjs) * parameterObj.block_length
            data_profiles_summary.append("\t".join([ \
                str(pair_idx), \
                str(len(profileObjs)), \
                str(bases), \
                "%.4f" % (sum(profileObjs)/bases).fixed, \
                "%.4f" % (sum(profileObjs)/bases).hetA, \
                "%.4f" % (sum(profileObjs)/bases).hetAB, \
                "%.4f" % (sum(profileObjs)/bases).hetB, \
                "%.4f" % (sum(profileObjs)/bases).invariant, \
                "%.4f" % (sum(profileObjs)/bases).multiallelic, \
                "%.4f" % (sum(profileObjs)/bases).missing 
                ]))
        fn_profiles_summary_tsv = parameterObj.profiles_summary_tsv_f
        with open(fn_profiles_summary_tsv, 'w') as fh_profiles_summary_tsv:
            fh_profiles_summary_tsv.write("\n".join(data_profiles_summary) + "\n")
        return fn_profiles_by_block_id_tsv, fn_profiles_summary_tsv
        
class WindowDataObj(object):
    __slots__ = ["windowObjs"]

    def __init__(self):
        self.windowObjs = OrderedDict()

    def __len__(self):
        return len(self.windowObjs)

    def add_windowObj(self, windowObj):
        self.windowObjs[windowObj.window_id] = windowObj

    def write_window_metrics_tsv(self, outprefix, pair_ids_by_pair_idx):
        # plot actual profiles of windows
        fn_window_metrics_tsv = "%s.window_metrics.tsv" % (outprefix)
        data_window_metrics_tsv = []
        data_window_metrics_tsv.append("#%s" % "\t".join(["window_id", "length", "span", "mean_block_density", "mean_sample_count", "mean_pair_count"]))

        fn_window_bed = "%s.window.bed" % (outprefix)
        data_window_bed = []
        data_window_bed.append("#%s" % "\t".join(["CHROM", "START", "END", "length", "count_samples", "count_pairs", "window_id"]))
        
        fn_window_profiles_tsv = "%s.window.profiles.tsv" % (outprefix)
        data_window_profiles_tsv = []
        data_window_profiles_tsv.append("#%s" % ("\t".join([ \
            "window_id", \
            "pair_ids", \
            "hetA", \
            "hetB", \
            "hetAB", \
            "fixed", \
            "multiallelic", \
            "missing", \
            "pi_A", \
            "pi_B", \
            "d_xy", \
            "f_st" \
            ])))
        with tqdm(total=len(self.windowObjs), desc="[%] ", ncols=200, unit_scale=True) as pbar:
            for window_id, windowObj in self.windowObjs.items():
                data_window_metrics_tsv.append(str(windowObj))
                for bedObj in windowObj.bedObjs:
                    data_window_bed.append("\t".join([str(bedObj), window_id]))
                for pair_idx, profile in windowObj.profileObj_by_pair_idx.items():
                    metrics = windowObj.metrics_by_pair_idx[pair_idx]
                    data_window_profiles_tsv.append("\t".join([ \
                        window_id, \
                        "_vs_".join(sorted([pair_ids_by_pair_idx[pair_idx][0], pair_ids_by_pair_idx[pair_idx][1]])), \
                        str(profile.get("hetA", 0)), \
                        str(profile.get("hetB", 0)), \
                        str(profile.get("hetAB", 0)), \
                        str(profile.get("fixed", 0)), \
                        str(profile.get("multiallelic", 0)), \
                        str(profile.get("missing", 0)), \
                        str(metrics.get("pi_A", 0)), \
                        str(metrics.get("pi_B", 0)), \
                        str(metrics.get("d_xy", 0)), \
                        str(metrics.get("f_st", 0)) \
                    ]))
                pbar.update()

        with open(fn_window_metrics_tsv, 'w') as fh_window_metrics_tsv:
            fh_window_metrics_tsv.write("\n".join(data_window_metrics_tsv))
        with open(fn_window_bed, 'w') as fh_window_bed:
            fh_window_bed.write("\n".join(data_window_bed))
        with open(fn_window_profiles_tsv, 'w') as fh_window_profiles_tsv:
            fh_window_profiles_tsv.write("\n".join(data_window_profiles_tsv))
        return fn_window_metrics_tsv, fn_window_bed, fn_window_profiles_tsv

class WindowObj(object):
    __slots__ = ["contig_id", "start", "end", "length", "span", "window_id", "profileObj_by_pair_idx", "metrics_by_pair_idx", "cov_by_sample_idx", "cov_by_pair_idx", "bedObjs"]

    def __init__(self, contig_id, blockObjs, block_length):
        self.contig_id = contig_id
        self.start = blockObjs[0].start
        self.end = blockObjs[-1].end
        self.length = block_length * len(blockObjs)
        self.span = blockObjs[-1].end - blockObjs[0].start
        self.window_id = "%s_%s_%s" % (contig_id, blockObjs[0].start, blockObjs[-1].end)
        self.profileObj_by_pair_idx = OrderedDict()
        self.metrics_by_pair_idx = OrderedDict()
        self.cov_by_sample_idx = {}
        self.cov_by_pair_idx = {}
        self.bedObjs = []
        self.populate(blockObjs)
        #self.compute_metrics()

    def __str__(self):
        return "%s" % "\t".join([\
                self.window_id, \
                str(self.length), \
                str(self.span), \
                "%.4f" % (self.length / self.span), \
                "%.4f" % (sum([cov for sample_idx, cov in self.cov_by_sample_idx.items()]) / (WINDOW_SIZE * SAMPLES_COUNT)), \
                "%.4f" % (sum([cov for pair_idx, cov in self.cov_by_pair_idx.items()]) / (WINDOW_SIZE * PAIRS_COUNT)) \
            ])

    def compute_metrics(self):
        for pair_idx, profile in self.profileObj_by_pair_idx.items():
            self.metrics_by_pair_idx[pair_idx] = calculate_metrics(profile, self.length)

    def populate(self, blockObjs):
        for blockObj in blockObjs:
            for pair_idx in blockObj.profileObj_by_pair_idx:
                self.profileObj_by_pair_idx[pair_idx] = self.profileObj_by_pair_idx.get(pair_idx, Counter()) + blockObj.profileObj_by_pair_idx[pair_idx]
                self.cov_by_pair_idx[pair_idx] = self.cov_by_pair_idx.get(pair_idx, 0) + 1
            for sample_idx in blockObj.sample_idxs:
                self.cov_by_sample_idx[sample_idx] = self.cov_by_sample_idx.get(sample_idx, 0) + 1
            for bedObj in blockObj.bedObjs:
                self.bedObjs.append(bedObj)

class BlockLengthException(Exception):
    pass

class PopsCountException(Exception):
    pass

if __name__ == "__main__":
    pass