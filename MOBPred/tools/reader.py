import os

class ReaderError(Exception):
    pass

class Mol2Format(object):
    """
    class used to read .mol2 file(s)
    """
    def open(self, filename, kwargs):
        import mol2
        return mol2.Reader(filename, **kwargs)

class PDBFormat(object):
    """
    class used to read .mol2 file(s)
    """
    def open(self, filename, kwargs):
        import pdb
        return pdb.Reader(filename, **kwargs)

known_formats={'.mol2': Mol2Format(), '.pdb': PDBFormat()}

class ReaderFormatError(ReaderError):

    def __init__(self, format, known_formats):
        self.format = format
        self.known_formats = known_formats

    def __str__(self):
        global known_formats
        return 'Unknown molecule file format "%s"\n Available formats are %s.\n' % (self.format, self.known_formats.keys())

def open(filename, **kwargs):

    format = os.path.splitext(filename)[1]
    if format not in known_formats:
        raise ReaderFormatError(format, known_formats)

    return known_formats[format].open(filename, kwargs)
