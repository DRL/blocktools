"""usage: blocktools fixcoordinates -y <FILE> -c <FILE> [-t <INT> -h]
    
    -h, --help
    -y, --yaml <FILE>                       YAML file with parameters
    -c, --coordinates <FILE>                Coordinates file 
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

def task_parse_coordinate_transformation_f(parameterObj):
    start = timer()
    print("[#] Parse coordinate transformation file ...")
    coordinateTransformObj = parse_coordinate_transformation_f(parameterObj)
    print("[+] Read file in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return coordinateTransformObj

def task_transform_coordinates(parameterObj, blockDataObj, coordinateTransformObj):
    start = timer()
    print("[#] Transforming coordinates ...")
    blockDataObj = transform_coordinates(parameterObj, blockDataObj, coordinateTransformObj)
    print("[+] Transformed coordinates in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return blockDataObj

def task_write_fixcoordinates_output(parameterObj, blockDataObj):
    start = timer()
    print("[#] Writing output ...")
    fn_block_bed, fn_block_bed_void = blockDataObj.write_block_bed(parameterObj)
    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_block_bed, timer() - start, memory_usage_psutil()))
    if fn_block_bed_void:
        print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_block_bed_void, timer() - start, memory_usage_psutil()))

def main():
    start_time = timer()
    args = docopt(__doc__)
    parameterObj = task_parse_parameters(args)
    #print(parameterObj.__dict__)
    blockDataObj = task_load_blockDataObj(parameterObj)
    coordinateTransformObj = task_parse_coordinate_transformation_f(parameterObj)
    blockDataObj = task_transform_coordinates(parameterObj, blockDataObj, coordinateTransformObj)
    task_write_fixcoordinates_output(parameterObj, blockDataObj)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass