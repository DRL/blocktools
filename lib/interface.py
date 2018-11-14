"""
blocktools

    Usage:
        ./blocktools <module> [<args>...]         

    Modules:
        blocks                 Makes blocks (block.bed) from multiinter BED file
        variants               Fetches and analyses variants based on VCF file and block.bed
        fixcoordinates         Relabels BED file of blocks according to coordinate file
        windows                Makes windows of blocks based on profile.tsv and    
        plot                   Makes plots


    Options:
        -h, --help                         Show this screen.
        -v, --version                      Show version.

    Help:
        https://blocktools.readme.io/ (TBD)
"""

import sys
from docopt import docopt
from timeit import default_timer as timer

def main():
    try:
        __version__ = '0.7.0'
        start_time = timer()
        args = docopt(__doc__, version=__version__, options_first=True)
        if args['<module>'] == 'blocks':
            import lib.commands.blocks as blocks
            blocks.main()
        elif args['<module>'] == 'variants':
            import lib.commands.variants as variants
            variants.main()
        elif args['<module>'] == 'windows':
            import lib.commands.windows as windows
            windows.main()
        elif args['<module>'] == 'plot':
            import lib.commands.plot as plot
            plot.main()
        elif args['<module>'] == 'fixcoordinates':
            import lib.commands.fixcoordinates as fixcoordinates
            fixcoordinates.main()
        else:
            sys.exit("%r is not a blocktools module. See 'blocktools -help'." % args['<module>'])
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)
