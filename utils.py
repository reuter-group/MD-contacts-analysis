import numpy as np
from numba import jit, prange
import time
import platform
import MDAnalysis as mda 

sel_glue_atoms = lambda group: ' name ' + ' '.join([atom for atom in group]) 
sel_glue_residues = lambda group: ' or '.join(['resname ' + r for r in group]) 


# @timing
def extract_coords(u, atom_group1, atom_group2, start=0, stop=None, step=1):
    if not isinstance(u, mda.Universe) or \
        not isinstance(atom_group1, mda.AtomGroup) or \
        not isinstance(atom_group2, mda.AtomGroup):
        raise TypeError

    M = len(atom_group1)
    N = len(atom_group2)
    U = len(u.trajectory[start:stop:step])
    coords1 = np.empty((U, M, 3), np.float32)
    coords2 = np.empty((U, N, 3), np.float32)

    for i, _ in enumerate(u.trajectory[start:stop:step]):
        coords1[i] = atom_group1.positions
        coords2[i] = atom_group2.positions

    return coords1, coords2



def read_resfile(resfile):
    """
    residue file format example:
    GLY48
    ILE49
    return: [('48', 'GLY'), ('49', 'ILE'), ...]
    """
    res_info = []   # [(resid, resname) ]
    with open(resfile) as f:
        for l in f.readlines():
            resid, resname = l.strip()[3:], l[:3]
            res_info.append((resid, resname))
    return  res_info


@jit(nopython=True, parallel=True, cache=True, nogil=True)
def filter_timeseries(arrs, n):
    if n <= 1:
        return None
    U, MN = arrs.shape
    for j in prange(MN):
        i = 0
        arr = arrs[:,j]
        while i < U:
            if arr[i] == 1:              
                if np.sum(arr[i: i+n]) == n: 
                    i += n
                elif i>=n and np.sum(arr[i-n: i]) == n:
                    i += 1
                else:
                    arr[i] = 0
                    i += 1
            else:
                i += 1
            

# Decorator for performace measurement
def timing(func):
    def clocked(*args):
        t0 = time.perf_counter()
        result = func(*args)
        elapsed = time.perf_counter() - t0
        name, module = func.__name__, func.__module__
        arg_str = ', '.join(str(arg)[:50] for arg in args) # print out only [:50]
        date = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        suffix = ' [at %s, %s]' % (platform.node(), date)

        REPORT = './performance_report.txt'  # report filename
        with open(REPORT, 'a') as f:
            f.write('[%0.4fs] %s.%s \nParameters: %s \n%s \n\n' %
                    (elapsed, module, name, arg_str, suffix))
        return result    
    return clocked
