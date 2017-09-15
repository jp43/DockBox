import numpy as np
import itertools as it

_known_entries = [ "ATOM  ", "HETATM", "ANISOU", "CRYST1",
    "COMPND", "MODEL", "ENDMDL", "TER", "HEADER", "TITLE", "REMARK",
    "CONECT"]

class PDBError(Exception):
    pass

class Reader(object):

    def __init__(self, filename, **kwargs):

        self.filename = filename
        self.file = open(filename, 'r')

    def read(self):

        struct = self.next()
        while struct is not None:
            yield struct
            struct = self.next()

    def next(self):

        struct = None
        for idx, line in enumerate(self.file):
            # initialize stucture
            if idx == 0:
                struct = {}
                struct['ATOM'] = []
            if line.startswith(('ATOM','HETATM')):
                line_s = [line[6:11], line[12:16], line[17:20], line[22:26], line[30:38], line[38:46], line[46:54]]
                struct['ATOM'].append(map(str.strip, line_s))
            elif line.startswith('END'):
                break
        return struct

    def _skip(self):
        try:
            self.file.next()
        except StopIteration:
            return None

        try:
            natoms = int(self.file.next())
        except StopIteration:
            raise PDBError("File ended unexpectedly when reading number of atoms.")

        for atom in it.izip(xrange(natoms), self.file):
            pass

        try:
            self.file.next()
        except StopIteration:
            raise PDBError("File ended unexpectedly when reading box line.")
        return None


    def readlines(self, *args):
        if len(args) == 0:
            configs = []
            config = self.next()
            while config is not None:
                configs.append(config)
                config = self.next()
        elif len(args) == 1:
            lines = args[0]
            if isinstance(lines, int):
                lines = [lines]
            else:
                lines = list(set(lines))
                lines.sort()
            lines = np.array(lines)
            lines = np.hstack((-1, lines))
            sklines = np.diff(lines) - 1
            configs = []
            for skline in sklines:
                for idx in xrange(skline):
                    self._skip()
                config = self.next()
                configs.append(config)
        else:
            raise PDBError("invalid number of arguments to readlines")

        return np.array(configs)

    def close(self):
        self.file.close()

    def __iter__(self):
        return self.read()

    readline = next

class Writer(object):
    pass
