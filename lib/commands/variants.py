"""usage: blocktools variants -y <FILE> [-t <INT> -h]
    
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

def task_load_blockDataObj(parameterObj):
    start = timer()
    print("[#] Loading blocks ...")
    blockDataObj = load_blockDataObj(parameterObj)
    print("[+] Read %s blocks in %.3fs (%.2fMB)" % (len(blockDataObj), timer() - start, memory_usage_psutil()))
    return blockDataObj

def task_fetch_variants(parameterObj, blockDataObj):
    start = timer()
    print("[#] Fetching variants ...")
    blockDataObj = parse_vcf(parameterObj, blockDataObj)
    print("[+] VCF parsed in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return blockDataObj

def task_analyse_variants(parameterObj, blockDataObj):
    start = timer()
    print("[#] Analysing variants ...")
    blockDataObj = analyse_variants(parameterObj, blockDataObj)
    print("[+] Analysed variants in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return blockDataObj

def task_write_variant_output(parameterObj, blockDataObj):
    start = timer()
    print("[#] Writing output ...")
    fn_profiles_by_block_id_tsv, profileObjs_by_pair_idx = blockDataObj.write_variant_blocks(parameterObj)
    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_profiles_by_block_id_tsv, timer() - start, memory_usage_psutil()))
    fn_profiles_summary_tsv = blockDataObj.write_variant_summary(parameterObj, profileObjs_by_pair_idx)
    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_profiles_summary_tsv, timer() - start, memory_usage_psutil()))

def main():
    start_time = timer()
    args = docopt(__doc__)
    parameterObj = task_parse_parameters(args)
    blockDataObj = task_load_blockDataObj(parameterObj)
    blockDataObj = task_fetch_variants(parameterObj, blockDataObj)
    blockDataObj = task_analyse_variants(parameterObj, blockDataObj)
    task_write_variant_output(parameterObj, blockDataObj)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass