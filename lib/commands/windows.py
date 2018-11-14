"""usage: blocktools windows -y <FILE> [-b <FILE> -t <INT> -w <INT> -l <INT> -h]
    
    -h, --help
    -y, --yaml <FILE>                       YAML file with parameters
    -b, --bed <FILE>                        BED file (if coordinates were fixed)
    -t, --threads <INT>                     Number of threads to use
    -w, --window_size <INT>                 Window size in blocks 
    -l, --overlap <INT>                     Window overlap in blocks
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

def task_load_profileObjs(parameterObj, blockDataObj):
    start = timer()
    print("[#] Loading profiles ...")
    blockDataObj = load_profileObjs(parameterObj, blockDataObj)
    print("[+] Read profiles in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return blockDataObj

def task_make_windows(parameterObj, blockDataObj):
    start = timer()
    print("[#] Making windows for %s contigs..." % (len(blockDataObj.blockObj_idxs)))
    windowDataObj = make_windows(parameterObj, blockDataObj)
    print("[+] Made %s windows in %.3fs (%.2fMB)" % (len(windowDataObj), timer() - start, memory_usage_psutil()))
    return windowDataObj

def task_analyse_windows(parameterObj, windowDataObj):
    start = timer()
    print("[#] Analysing windows ...")
    windowDataObj = analyse_windows(parameterObj, windowDataObj)
    print("[+] Analysed windows in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return windowDataObj

def task_write_window_output(parameterObj, windowDataObj):
    start = timer()
    print("[#] Writing output ...")
    fn_window_metrics_tsv, fn_window_bed, fn_window_profiles_tsv, fn_window_sfs_tally = windowDataObj.write_window_output(parameterObj)
    print("[+] Wrote \n\t'%s' and \n\t'%s' and \n\t'%s' and \n\t'%s' in %.3fs (%.2fMB)" % (fn_window_metrics_tsv, fn_window_bed, fn_window_profiles_tsv, fn_window_sfs_tally, timer() - start, memory_usage_psutil()))

def main():
    start_time = timer()
    args = docopt(__doc__)
    parameterObj = task_parse_parameters(args)
    blockDataObj = task_load_blockDataObj(parameterObj)
    blockDataObj = task_load_profileObjs(parameterObj, blockDataObj)
    windowDataObj = task_make_windows(parameterObj, blockDataObj)
    windowDataObj = task_analyse_windows(parameterObj, windowDataObj)
    task_write_window_output(parameterObj, windowDataObj)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass