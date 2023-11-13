#!/usr/bin/env python
import sys, os
import time
import math 
import json 
import psutil

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from numba import jit, prange
from hydrop_candi import ATOMS_DICT
from utils import *
 

############# Analysis Specific Configurations #############
DIS_CUTOFF = 3 # distance cutoff 
PRO_SEGID = 'PROA' # segid for protein residues
MEMB_SEGID = 'MEMB' # segid for memb lipids
############################################################


## Code can be executed in parallel
@jit(nopython=True, parallel=True, cache=True, nogil=True)
def cal_dis_mat(fnumber, M, N, atoms1, atoms2):
    mat = np.zeros((fnumber, M, N), dtype=np.bool_)
    for f in prange(fnumber):
        for i in prange(M):
            for j in prange(N):        
                # dis = np.sqrt(np.sum( (atoms1[f][i]- atoms2[f][j])**2))
                dis = ((atoms1[f][i][0] - atoms2[f][j][0])**2 + 
                        (atoms1[f][i][1] - atoms2[f][j][1])**2 +
                        (atoms1[f][i][2] - atoms2[f][j][2])**2 )**.5
                mat[f][i][j] = dis < DIS_CUTOFF
    return mat


def sub_ana(group1, group2, u, lifetime=None, details=False):
    ## Extract coordinates of all candidate atoms from each frame
    t3 = time.perf_counter()
    coor1, coor2 = extract_coords(u, group1, group2, fstart, fstop, fstep) 

    # # Calculate the distance matrix of all candidates and keep those that meet DIS_CUTOFF
    t4 = time.perf_counter()
    N_sub = len(group2)
    times_matrix = cal_dis_mat(U, M, N_sub, coor1, coor2) 
    
    # Remove those that occurred less than <lifetime> consecutive frames 
    t5 = time.perf_counter()
    if lifetime is not None:
        filter_timeseries(times_matrix.reshape(U, -1), lifetime) # U * (M*N)

    t6 = time.perf_counter()
    # Count the contacts for each amino acid residue
    offset = 0
    for residue_atoms in atom_list1:
        key = "{}-{}".format(residue_atoms[0].resid, residue_atoms[0].resname.upper())
        res_mat = times_matrix[:, offset: offset+len(residue_atoms),:] 
        if details:
            for i, frame in enumerate(res_mat):
                ixs = np.nonzero(frame)
                for j, ix in enumerate(ixs[0]):
                    atomP = residue_atoms[ix]
                    atomM = group2[ixs[1][j]]
                    contact_dict[key].append((
                        i+1,
                        atomP.name, atomP.id, 
                        atomM.resname, atomM.resid, atomM.name, atomM.id))
        else:
            contact_dict[key] += np.sum(res_mat)   # only number of contacts 

        offset += len(residue_atoms)

    return (t4-t3, t5-t4, t6-t5)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: please submit a configuration file. Usage:\n ./hydrophobic_ana_v2.py CONFIG_FILE")
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
    LIPIDS = Parameters['hydrophobic_lipids']

    t0 = time.perf_counter()
    # Select interested candidate atoms in protein
    atom_group1 = mda.AtomGroup([],u)
    atom_list1 = []
    contact_dict = {}
    print("Selecting protien atoms from segid {}...".format(PRO_SEGID))
    for rid, rname in res_info:
        if rname.upper() in ATOMS_DICT:
            sel_prot = '(segid {} and resid {} and ({}))'.format(PRO_SEGID, rid, 
                sel_glue_atoms(ATOMS_DICT[rname.upper()].split()))
            temp_group = u.select_atoms(sel_prot)
            atom_group1 += temp_group
            atom_list1.append(temp_group)
            key = "{}-{}".format(rid, rname.upper())
            contact_dict[key] = 0 if details_file == "" else []

    t1 = time.perf_counter()
    # Select interested candidate atoms in lipids
    print("Selecting membrane atoms from segid {}...".format(MEMB_SEGID))
    sel_memb = 'segid {} and '.format(MEMB_SEGID)
    for lipid in LIPIDS:
        if lipid in ATOMS_DICT:
            sel_memb += ' ( resname {} and ({}) ) '.format(
                lipid, 
                sel_glue_atoms(set(ATOMS_DICT[lipid].split()) - set(ATOMS_DICT[lipid + '_HEAD'].split()))
            )
            if lipid != LIPIDS[-1]:
                sel_memb += 'or'
    atom_group2 = u.select_atoms(sel_memb) 

    t2 = time.perf_counter()
    '''
    This selection overhead scales linealy with the number of interested amino acids and lipids
    but independently from the number of frames
    '''
    print("*Selection time usage: protein {:.2f}s, lipids {:.2f}s".format(t1-t0, t2-t1))

    # Print out input size information
    M = len(atom_group1)
    N = len(atom_group2)
    U = len(u.trajectory[fstart:fstop:fstep])
    print("\nTotal {} frames, selected {}: {} atoms, {}: {} atoms\n".format(U, PRO_SEGID , M , MEMB_SEGID, N))

    ## Split lipid groups into smaller subgroups which can fit the memory limitation
    # note: splitting by frames will lead to dealing with connection issue of subtrajectories
    # when applying lifetime filtering, which will be tricky
    mem_avail = psutil.virtual_memory().available  # system's available memory for starting new applications
    mem_bytes = min(mem_avail, Parameters['mem_limit'] * 1024**3)
    coor1_mem_size = M * U * 3 * 4 # bytes
    s_size = math.floor((mem_bytes - coor1_mem_size) / U / M )  # total number of atoms in each subgroup
    s_num = math.ceil( N / s_size ) # total number of subgroups 
    subgroups = [atom_group2[i*s_size : (i+1)*s_size] for i in range(s_num)]  

    ## Run the calculation subgroup by subgroup
    acct1, acct2, acct3 = 0, 0, 0
    for i, subgroup in enumerate(subgroups):
        print("[{:.2f}%] Running...".format(i /len(subgroups) * 100 ))
        t10, t20, t30 = sub_ana(atom_group1, subgroup, u, lifetime=lifetime, details= details_file!="")
        acct1 += t10
        acct2 += t20
        acct3 += t30 
     
    print("Time usage: \n extract coordinates: {:.2f}s \n calculate distance matrix: {:.2f}s \
            \n lifetime cutoff: {:.2f}s".format(acct1, acct2, acct3)) 
    
    if details_file == "":
        for k,v in contact_dict.items():
            print("{} average contacts per frame: {:.2f}".format(k, v/U))
    else:
        for k,v in contact_dict.items():
            print("{} average contacts per frame: {:.2f}".format(k, len(v)/U))
            if len(v) > 0:
                rid, rname = k.split('-')
                with open(details_file+'_{}_hydrophobic_contacts.txt'.format(k), 'w') as f:
                    sorted_v = sorted(v, key=lambda i: (i[0], i[2])) # sort by frame index
                    for c in sorted_v :
                        f.write("[{}] {} {} : {} {} -- {} {} : {} {}\n".format(
                            c[0], rname, rid, c[1], c[2], c[3], c[4], c[5], c[6]))
                                            
