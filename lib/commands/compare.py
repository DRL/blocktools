"""usage: blocktools compare -1 <FILE> -2 <FILE> [--p1 <STR> --p2 <STR> -h]
    
    -h, --help
    -1, --yaml1 <FILE>                       YAML file 1 with parameters
    --p1 <STR>                               Alternative outprefix 1
    -2, --yaml2 <FILE>                       YAML file 2 with parameters
    --p2 <STR>                               Alternative outprefix 2
"""

from docopt import docopt
from lib.src.code import *
from timeit import default_timer as timer

def task_parse_parameters(args):
    start = timer()
    print("[#] Parsing parameters ...")
    parameterObj = parse_parameters(args)
    print("[+] Read parameters in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return parameterObj

def task_parse_genome_f(parameterObj):
    start = timer()
    print("[#] Parsing genome file ...")
    sequence_OrdDict = parse_genome_f(parameterObj)
    print("[+] Read %s sequences in %.3fs (%.2fMB)" % (len(sequence_OrdDict), timer() - start, memory_usage_psutil()))
    return sequence_OrdDict

def task_plot_pairs_comparison(parameterObj_1, parameterObj_2):
    start = timer()
    print("[#] Parsing %s ..." % parameterObj_1.variant_pairs_tsv_f)
    pairs_fst_comparison_fn = plot_pairs_comparison(parameterObj_1, parameterObj_2)
    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (pairs_fst_comparison_fn, timer() - start, memory_usage_psutil()))

def main():
    start_time = timer()
    args = docopt(__doc__)
    args_1 = {'--yaml' : args['--yaml1'], '--outprefix': args['--p1']}
    args_2 = {'--yaml' : args['--yaml2'], '--outprefix': args['--p2']}
    parameterObj_1 = task_parse_parameters(args_1)
    parameterObj_2 = task_parse_parameters(args_2)
    task_plot_pairs_comparison(parameterObj_1, parameterObj_2)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass