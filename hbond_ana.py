#!/usr/bin/env python
import numpy as np
import sys
import MDAnalysis as mda
import MDAnalysis.analysis.hbonds
from hbond_candi import *
from utils import filter_timeseries, read_resfile
import time
import json 

########## Modify these configurations if necessary ########
DIS_CUTOFF = 2.4   # distance cutoff between H and A
ANGLE_CUTOFF = 130.0, # angle cutoff of D--H--A
PRO_SEGID = 'PROA' # segid for protein residues
MEMB_SEGID = 'MEMB' # segid for memb lipids
##### NOTE for Lipid Types #####
# the lists of donors/acceptors of all lipids will be merged
# and used as candidate donors/acceptors to run the analysis
############################################################


def _set_lipid_atoms(lipids, htype='donor'):
    '''
    Merge the donors/acceptors of all the given lipids
    '''
    atoms = []
    if htype == 'donor':
        for l in lipids:
            if l in ATOMS_DICT:            
                atoms += ATOMS_DICT[l]['donor'].split()
    elif htype == 'acceptor':
        for l in lipids:
            if l in ATOMS_DICT:            
                atoms += ATOMS_DICT[l]['acceptor'].split()
    return atoms



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: please submit a configuration file. Usage:\n ./hbond_ana.py CONFIG_FILE")
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
    LIPIDS = Parameters['hbond_lipids']

    rids = [r[0] for r in res_info]
    t0 = time.perf_counter()

    h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(
        u, 
        selection1 = 'segid {} and resid {}'.format(PRO_SEGID, ' '.join(rids) ) , 
        selection2 = 'segid {} and resname {}'.format(MEMB_SEGID, ' '.join(LIPIDS) ) , 
        selection1_type='both', # consider sel1 as both donor and acceptor
        update_selection1=False,  # candidates no need to change in each frame
        update_selection2=False, 
        distance = DIS_CUTOFF,   
        angle = ANGLE_CUTOFF, 
        filter_first = False, # a coarse distance fileter; not necessary 
        
        ### further selections for donors/acceptors performed on sel1/sel2
        donors= _set_lipid_atoms(LIPIDS, 'donor'), 
        acceptors= _set_lipid_atoms(LIPIDS, 'acceptor'),
    )
    t1 = time.perf_counter()  # HydrogenBondAnalysis setting-up time

    h.run(fstart, fstop, fstep)  # run the analysis
    t2 = time.perf_counter()
    
    hbonds = h.count_by_type()  # shorter list including each type of hbond
    recs = h.timesteps_by_type() # longer list plus timestep information
    # print(hbonds)
    # print(recs)
    ### NOTE: type of recs 
    #  dtype = [
    #     ('donor_index', int),
    #     ('acceptor_index', int), ('donor_resnm', 'U4'), ('donor_resid', int),
    #     ('donor_heavy_atom', 'U4'), ('donor_atom', 'U4'),('acceptor_resnm', 'U4'),
    #     ('acceptor_resid', int), ('acceptor_atom', 'U4'), ('time', float)]

    ### Post-processing for 
    # 1. perfrom lifetime cutoff selection
    # 2. customize the output format
    
    U = len(u.trajectory[fstart:fstop:fstep])
    frame_flags = np.zeros((len(hbonds), U ))
    for i, hbond in enumerate(hbonds): 
        for r in recs:
            if (r['donor_index'], r['acceptor_index']) == \
                  (hbond['donor_index'], hbond['acceptor_index']):
                # convert timestep to frame index
                frame_flags[i][h.timesteps.index(r['time'])] = 1  

    # perform lifetime cutoff
    t3 = time.perf_counter()
    if lifetime is not None:
        filter_timeseries(frame_flags.transpose() , lifetime)
    t4 = time.perf_counter()

    print("\nTotal {} frames".format(U))
    print("Time usage: \n - initialization(including atom seletion):{:.2f} \
           \n - finding out candidate hbonds: {:.2f}s \
           \n - lifetime cutoff: {:.2f}s".format(t1-t0, t2-t1, t4-t3)) 

    # store frame indices where hbond occurs for calculating occupancy
    new_hbonds = [] 
    for i, hbond in enumerate(hbonds):
        frames = frame_flags[i].nonzero()[0] 
        new_hbonds.append(list(hbond)[:9] + [frames])

    ### Output the results 
    for rid, rname in res_info :
        head_idx, gly_idx, pho_idx = [], [], []
        for hbond in new_hbonds:
            if hbond[2] == rname and hbond[3] == int(rid): # amino acid is donor
                lipid_name, lipid_atom = hbond[6], hbond[8]
            elif hbond[6] == rname and hbond[7] == int(rid): # amino acid is acceptor
                lipid_name, lipid_atom = hbond[2], hbond[4]
            else:
                continue
    
            if lipid_atom in ATOMS_DICT[lipid_name]['pho'].split():
                pho_idx += list(hbond[-1])
            elif lipid_atom in ATOMS_DICT[lipid_name]['gly'].split():
                gly_idx += list(hbond[-1])
            elif lipid_atom in ATOMS_DICT[lipid_name]['head'].split():
                head_idx += list(hbond[-1])           
            else:
                print("WARNING: found contact from parts other than PHO, GLY and HEAD")
                print(str(hbond))
    
        gly_occ = len(set(gly_idx)) / U * 100
        head_occ = len(set(head_idx)) / U * 100
        pho_occ = len(set(pho_idx)) / U * 100
        
        print("{}-{} Hbond occupancy: PHOSPHATE {:.2f}%   GLYCEROL {:.2f}%   HEADGROUP {:.2f}% "\
            .format(rname, rid, pho_occ, gly_occ, head_occ))

    # output details to file
    if details_file != "":
        with open(details_file+'_hbonds.txt', mode='w') as f:
            f.write(','.join([
                        'RES-ID',
                        'donor_index', 'acceptor_index', 'donor_resnm', 
                        'donor_resid', 'donor_heavy_atom', 'donor_atom', 
                        'acceptor_resnm','acceptor_resid', 'acceptor_atom', 
                        'frame_index'
                        ]) 
                    + '\n')

            for rid, rname in res_info:
                for hbond in new_hbonds:
                    if hbond[2] == rname and hbond[3] == int(rid) or \
                       hbond[6] == rname and hbond[7] == int(rid): 
                        for f_number in hbond[-1]:
                            strs = [str(field) for field in hbond[:9]]
                            strs.insert(0, rname+'-'+rid)
                            strs.append(str(f_number))
                            f.write(' '.join(strs) + '\n')
