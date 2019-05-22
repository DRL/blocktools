"""usage: blocktools unpickle -i <FILE> [-h]
    
    -h, --help
    -i, --pickle <FILE>                     Pickle file
"""

from docopt import docopt
import pickle

def main():
    args = docopt(__doc__)
    print(pickle.load(open(args['--pickle'], "rb")))

if __name__ == "__main__":
    pass
