import os

class WriterError(Exception):
    pass

class Mol2Format(object):
    """
    class used to write files with a single column
    """
    def open(self, *kwargs):
        import mol2
        return mol2.Writer()


known_formats={'.mol2': Mol2Format()}

class WriterFormatError(WriterError):

    def __init__(self, format, known_formats):
        self.format = format
        self.known_formats = known_formats

    def __str__(self):
        global known_formats
        return 'Unknown molecule file format "%s"\n Available formats are %s.' % (self.format, self.known_formats.keys())


def open(format, **kwargs):
    """
    writer.open('.eps') will create a sl.writer instance
    """

    if format not in known_formats:
        raise WriterFormatError(format, known_formats)

    return known_formats[format].open(kwargs)
