#!/usr/bin/env python
import numpy as np
import MDAnalysis as mda 
from utils import *
import sys
import time
import json


########## Modify these configurations if necessary #########
DIS_CUTOFF = 7 # distance cutoff between each C atom on the ring and N atom in the lipid
DIS_DIFF_CUTOFF = 1.5  # distances should not differ by more than 1.5 A
PRO_SEGID = 'PROA' # segid for protein residues
MEMB_SEGID = 'MEMB' # segid for memb lipids
AROMATIC = { # the aromatic ring atoms for candidate residues, from CHARMM script 'timeseries.inp'
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TRP': ['CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'NE1'],
}
############################################################


## Code can be executed in parallel
@jit(nopython=True, parallel=True, cache=True, nogil=True)
def cal_dis_mat(fnumber, M, N, atoms1, atoms2, Rix):
    '''
    fnumber: frame number, 
    atoms1.shape: (fnumber, M, 3)
    atoms2.shape: (fnumber, N, 3)
    Return: np.ndarray, (fnumber, M, len(Rix))
    ''' 
    mat = np.zeros((fnumber, M, len(Rix)), dtype=np.bool_)
    for f in prange(fnumber):
        for i in prange(M):
            dis_row = np.zeros((N,), dtype=np.float32)
            for j in prange(N):        
                # dis_row[j] = np.sqrt(np.sum( (atoms1[f][i]- atoms2[f][j])**2))
                dis_row[j] = ((atoms1[f][i][0] - atoms2[f][j][0])**2 + 
                                (atoms1[f][i][1] - atoms2[f][j][1])**2 +
                                (atoms1[f][i][2] - atoms2[f][j][2])**2 )**.5
                
            for r in prange(len(Rix)):                
                ring = dis_row[Rix[r][0]:Rix[r][1]]
                mat[f][i][r] = (ring < DIS_CUTOFF).all() and max(ring) - min(ring) < DIS_DIFF_CUTOFF
    return mat


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: please submit a configuration file. Usage:\n ./cation_pi_ana.py CONFIG_FILE")
        sys.exit(0)
        
    config_file = sys.argv[1]
    with open(config_file) as f:
        Parameters = json.load(f)
    
    # Set up the parameters
    u = mda.Universe(Parameters['PSF_FILE'] , Parameters['DCD_FILE'])
    res_info = read_resfile(Parameters['RESIDUE_FILE'])
    lifetime = Parameters.get('lifetime')
    details_file = Parameters.get('output_file_prefix').strip() 
    fstart, fstop, fstep = Parameters['start'], Parameters['stop'], Parameters['step']
    LIPIDS = Parameters['cation_pi_lipids']
    
    # Select lipid candidate atoms
    sel_N = 'segid {} and ({}) and name N'.format(MEMB_SEGID, sel_glue_residues(LIPIDS))
    atomsN = u.select_atoms(sel_N)  
    
    # Select residue candidate atoms 
    atomsC = mda.AtomGroup([],u)
    atomsC_ix = []  
    index0 = 0
    for rid, rname in res_info:
        Cs = AROMATIC.get(rname.upper())
        if Cs is None:
            print("WARNING: {}-{} is not considered as candidate.".format(rname, rid))
            continue
        sel_C = 'segid {} and resid {} and ({})'.format(PRO_SEGID, rid, sel_glue_atoms(Cs))
        temp_group = u.select_atoms(sel_C)  # aromatic ring atoms   
        atomsC += temp_group
        atomsC_ix.append([index0, index0+len(temp_group)]) # indices of residue's start and end 
        index0 += len(temp_group)

    # Print out input size information
    C = len(atomsC)
    N = len(atomsN)
    U = len(u.trajectory[fstart:fstop:fstep])
    print("\nTotal {} frames, selected {}: {} atoms, {}: {} atoms".format(U, PRO_SEGID , C , MEMB_SEGID, N))

    # Extract coordinates of all candidate atoms from each frame
    t0 = time.perf_counter()
    coorC, coorN = extract_coords(u, atomsC, atomsN, fstart, fstop, fstep)

    # Calculate the distances and only keep those that meet the cutoffs
    t1 = time.perf_counter()    
    times_matrix = cal_dis_mat(U, N, C, coorN, coorC, np.array(atomsC_ix))  
    ## shape: U * N * R (R = len(res_info))

    # Remove the contacts that not meet lifetime cutoff   
    t2 = time.perf_counter()
    if lifetime is not None:
        filter_timeseries(times_matrix.reshape(U, -1), lifetime) 
    
    t3 = time.perf_counter()
    print("Time usage: \n extract coordinates: {:.2f}s \n calculate distance matrix: {:.2f}s \
           \n lifetime cutoff: {:.2f}s\n".format(t1-t0, t2-t1, t3-t2)) 
    
    ## Output summary (and details to file)
    for i in range(len(res_info)):
        rid, rname = res_info[i]  
        res_mat = times_matrix[:,:,i].reshape(U,N)        
        occ = np.count_nonzero(np.sum(res_mat, axis=1) ) / U
        print("{}-{}: occupancy {:.2f}%".format(rname, rid, occ*100 ))

        if details_file != "":
            if occ > 0:
                with open(details_file+'_{}{}_cation_pi_contacts.txt'.format(rname, rid), 'w') as f:                
                    ixs = np.nonzero(res_mat)
                    for j, ix in enumerate(ixs[0]):
                        atomN = atomsN[ixs[1][j]]
                        f.write("[{}] {} {}  {} {}\n".format(ix+1, atomN.resname, atomN.resid, atomN.name, atomN.id)) 
    

## Note: results comparison for TYR44, last_50ns (no lifetime cutoff)
#contacted lipid_id      POPC-98  POPC-93   POPC-44   Total contacts/frame
#results of this script    698      187      24        909/2501 = 0.363
#results of old script     701      190      missing   891/2501 = 0.356