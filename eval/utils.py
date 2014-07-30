import sys
import pdb
import scipy as sp

def read_fasta(filename):
    """Reads a sequence in FASTA format from filename and stores it
       in a dictionary - one entry per contig"""

    genome = dict()

    for line in open(filename, 'r'):
        if line[0] == '>':
            if not is_first:
                genome[chr] = ''.join(sequence)
            is_first = False
            chr = line.strip()[1:]
            print 'Reading %s' % chr
            sequence = []
        else:
            sequence.append(line.strip())

    if len(sequence) > 0:
        genome[chr] = ''.join(sequence)

    return genome

def dict2hdf5(filename, in_dict):
    """Stores dictionary in hdf5 format. Works recursively."""

    import h5py

    if isinstance(filename, str): 
        ### open outfile
        outfile = h5py.File(filename, 'w')
        for key in in_dict:
            if isinstance(in_dict[key], dict):
                grp = outfile.create_group(name=key)
                dict2hdf5(grp, in_dict[key])
            else:
                outfile.create_dataset(name=key, data=in_dict[key], chunks=True, compression='gzip')
        outfile.close()
    else:
        for key in in_dict:
            if isinstance(in_dict[key], dict):
                grp = filename.create_group(name=key)
                dict2hdf5(grp, in_dict[key])
            else:
                filename.create_dataset(name=key, data=in_dict[key], chunks=True, compression='gzip')

def hdf52dict(filename):
    
    import h5py

    rt_dict = dict()

    IN = h5py.File(filename, 'r')
    for key in IN:
        rt_dict[key] = IN[key][...]
    IN.close()

    return rt_dict


def unique_rows(x):
    """This function takes a 2D scipy array x and makes it unique by rows."""
    
    y = sp.ascontiguousarray(x).view(sp.dtype((sp.void, x.dtype.itemsize * x.shape[1])))
    _, idx = sp.unique(y, return_index = True)

    return x[idx]


def get_mm(sl):
    """Gets an alignment line in SAM format and returns the number of mismatches"""

    for opt in sl[11:]:
        if opt[:3] == 'NM:':
            return int(opt[5:])
    return -1
