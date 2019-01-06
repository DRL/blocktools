"""usage: blocktools blocks -y <FILE> [-t <INT> -h]
    
    -h, --help
    -y, --yaml <FILE>                       YAML file with parameters
    -t, --threads <INT>                     Number of threads to use

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

def task_parse_multibed_f(parameterObj, sequence_OrdDict):
    start = timer()
    print("[#] Loading multiBED features ...")
    regionBatchObjs, coverageObj = parse_multibed_f(parameterObj, sequence_OrdDict)
    print("[+] Intervals read in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return regionBatchObjs, coverageObj

def task_make_blocks(parameterObj, regionBatchObjs):
    start = timer()
    print("[#] Cutting blocks using algorithm %s ..." % parameterObj.algorithm)
    blockDataObj = make_blocks(parameterObj, regionBatchObjs)
    print("[+] Made %s blocks in %.3fs (%.2fMB)" % (len(blockDataObj), timer() - start, memory_usage_psutil()))
    return blockDataObj

def task_write_block_output(parameterObj, blockDataObj, coverageObj):
    start = timer()
    print("[#] Writing output ...")
    fn_block_bed, fn_block_void_bed = blockDataObj.write_block_bed(parameterObj)
    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_block_bed, timer() - start, memory_usage_psutil()))
    if fn_block_void_bed:
        print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_block_void_bed, timer() - start, memory_usage_psutil()))
    start = timer()
    fn_block_tsv = blockDataObj.write_block_pairs(parameterObj)
    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_block_tsv, timer() - start, memory_usage_psutil()))
    start = timer()
    fn_span_tsv = blockDataObj.write_block_summary(parameterObj)
    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_span_tsv, timer() - start, memory_usage_psutil()))
    start = timer()
    print("[#] Writing output ...")
    print(vars(coverageObj))
    write_yaml(vars(coverageObj), parameterObj.bed_coverage_f)
    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (parameterObj.bed_coverage_f, timer() - start, memory_usage_psutil()))

def main():
    start_time = timer()
    args = docopt(__doc__)
    #print(args)
    parameterObj = task_parse_parameters(args)
    #print(parameterObj.__dict__)
    sequence_OrdDict = task_parse_genome_f(parameterObj)
    regionBatchObjs, coverageObj = task_parse_multibed_f(parameterObj, sequence_OrdDict)
    blockDataObj = task_make_blocks(parameterObj, regionBatchObjs)
    task_write_block_output(parameterObj, blockDataObj, coverageObj)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass