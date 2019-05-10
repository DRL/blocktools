# -*- coding: utf-8 -*-

from os import makedirs
from os.path import join, isdir

def create(out_f, path, header, lines):
    if path and not isdir(path):
        makedirs(path, exist_ok=True)
    if out_f:
        if lines:
            out_f = join(path, out_f)
            with open(out_f, 'w') as out_fh:
                if header:
                    out_fh.write("%s\n" % header)
                out_fh.write("%s\n" % "\n".join(lines))
        print("[>]\tCreated: %r" % out_f)
    else:
        print("[>]\tCreated: %r" % path)