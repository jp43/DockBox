import numpy as np

def center_of_geometry(coords):
    cog = np.array([0.0, 0.0, 0.0])
    coords = np.array(coords)
    num_atom = coords.shape[0]
    for idx in xrange(num_atom):
        cog = cog + coords[idx]
    cog = cog / num_atom
    return cog

