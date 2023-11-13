'''
Customized donors/acceptors for lipids
'''
ATOMS_DICT = { 
   
   ### *PE lipids
   'POPE': {
      'donor' : 'N HN1 HN2 HN3',
      'acceptor' : 'O13 O14 O12 O11 O21 O22 O31 O32',
      'pho': 'O13 O14 O12 O11', # phosphate
      'gly': 'O21 O22 O31 O32', # glycerol
      'head': 'N HN1 HN2 HN3',  # head group
   },

   ### *PC lipids
   'POPC': {
      'donor': '', # no candidate donor for *PC
      'acceptor': 'O13 O14 O12 O11 O21 O22 O31 O32',
      'pho': 'O13 O14 O12 O11',
      'gly': 'O21 O22 O31 O32',
      'head': '', # no head for *PC
   },

   ### POPI lipids
   'POPI': {
      'donor' : 'O2 O3 O4 O5 O6',
      'acceptor' : 'O2 O3 O4 O5 O6 O13 O14 O12 O11 O21 O22 O31 O32',
      'pho': 'O13 O14 O12 O11', # phosphate
      'gly': 'O21 O22 O31 O32', # glycerol
      'head': 'O2 O3 O4 O5 O6',  # head group
   },

   'POPS': {
      'donor' : 'N',
      'acceptor' : 'O13 O14 O12 O11 O21 O22 O31 O32 O13A O13B N',
      'pho': 'O13 O14 O12 O11', # phosphate
      'gly': 'O21 O22 O31 O32', # glycerol
      'head': 'N O13A O13B',  # head group
   },
   
   'CER180': {
      'donor' : 'NF O3 O1',
      'acceptor' : 'OF NF O3 O1',
      'pho': '', # no phosphate for CER180
      'gly': '', # no glycerol for CER180
      'head': 'NF O3 O1 OF',  # head group
   },
   
   'CER1': {  # duplicate of 'CER180'
      'donor' : 'NF O3 O1',
      'acceptor' : 'OF NF O3 O1',
      'pho': '', # no phosphate for CER180
      'gly': '', # no glycerol for CER180
      'head': 'NF O3 O1 OF',  # head group
   },
   
   'CHL1': {
      'donor' : 'O3',
      'acceptor' : 'O3',
      'pho': '', # no phosphate for CHL1
      'gly': '', # no glycerol for CHL1
      'head': 'O3',  # head group
   },

   'PSM': {
      'donor' : 'NF O3',
      'acceptor' : 'O11 O12 O13 O14 OF',
      'pho': 'O11 O12 O13 O14',
      'gly': 'OF NF O3', 
      'head': '',  # no candidate
   },

   'SAPI24':{  # 
      'donor' : 'O2 O3 O6 OP42',
      'acceptor' : 'O21 O22 O31 O32 O11 O12 O13 O14 O5 OP54 OP53 OP52 O4 OP43 OP44',
      'pho': 'O11 O12 O13 O14',
      'gly': 'O21 O22 O31 O32', 
      'head': 'O2 O3 O6 O4 OP42 OP43 OP44 O5 OP52 OP53 OP54'
   },

   # In case long lipid names (> 4 letters) can't be recognised:
   # we just duplicate the entry for 'SAPI24' and rename it with the first 4 letters
   'SAPI':{  # 
      'donor' : 'O2 O3 O6 OP42',
      'acceptor' : 'O21 O22 O31 O32 O11 O12 O13 O14 O5 OP54 OP53 OP52 O4 OP43 OP44',
      'pho': 'O11 O12 O13 O14',
      'gly': 'O21 O22 O31 O32', 
      'head': 'O2 O3 O6 O4 OP42 OP43 OP44 O5 OP52 OP53 OP54'
   },

}

'''
.. Default donors/acceptors used in MDAnalysis.analysis.hbonds.HydrogenBondAnalysis 
    (CHARMM27 force field)

   Donors(15): 'N', 'OH2', 'OW', 'NE', 'NH1', 'NH2', 'ND2', 'SG', 'NE2', 'ND1', 'NZ', 'OG', 'OG1', 'NE1', 'OH',
   Acceptors(17): 'O', 'OC1', 'OC2', 'OH2', 'OW', 'OD1', 'OD2', 'SG', 'OE1', 'OE1', 'OE2', 'ND1', 'NE2', 'SD', 'OG', 'OG1', 'OH',

.. : More details 
   =========== ==============  =========== ====================================
   group       donor           acceptor    comments
   =========== ==============  =========== ====================================
   main chain  N               O, OC1, OC2 OC1, OC2 from amber99sb-ildn
                                           (Gromacs)
   water       OH2, OW         OH2, OW     SPC, TIP3P, TIP4P (CHARMM27,Gromacs)

   ARG         NE, NH1, NH2
   ASN         ND2             OD1
   ASP                         OD1, OD2
   CYS         SG
   CYH                         SG          possible false positives for CYS
   GLN         NE2             OE1
   GLU                         OE1, OE2
   HIS         ND1, NE2        ND1, NE2    presence of H determines if donor
   HSD         ND1             NE2
   HSE         NE2             ND1
   HSP         ND1, NE2
   LYS         NZ
   MET                         SD          see e.g. [Gregoret1991]_
   SER         OG              OG
   THR         OG1             OG1
   TRP         NE1
   TYR         OH              OH
   =========== ==============  =========== ====================================
'''
