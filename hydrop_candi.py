# Inventory of the candidate atoms for the hydrophobic interaction analysis.  
# The atom names correspond to the nomenclature of the CHARMM36 force field

ATOMS_DICT = { 
    'ALA': 'CA HA CB HB1 HB2 HB3',
    'ARG': 'CA HA CB HB1 HB2 CG HG1 HG2',
    'ASN' : 'CA HA CB HB1 HB2',
    'ASP' : 'CA HA CB HB1 HB2',
    'CYS' : 'CA HA CB HB1 HB2',
    'GLN' : 'CA HA CB HB1 HB2 CG HG1 HG2',
    'GLU' : 'CA HA CB HB1 HB2 CG HG1 HG2',
    'GLY' : 'CA HA1 HA2',
    'HSD' : 'CA HA CB HB1 HB2',
    'HSE' : 'CA HA CB HB1 HB2',
    'HSP' : 'CA HA CB HB1 HB2 CG' ,
    'ILE' : 'CA HA CB HB CG1 HG11 HG12 CG2 HG21 HG22 HG23 CD HD1 HD2 HD3',
    'LEU' : 'CA HA CB HB1 HB2 CG HG CD1 HD11 HD12 HD13 CD2 HD21 HD22 HD23',
    'LYS' : 'CA HA CB HB1 HB2 CG HG1 HG2 CD HD1 HD2 CE HE1 HE2',
    'MET' : 'CA HA CB HB1 HB2 CG HG1 HG2 CE HE1 HE2 HE3',
    'PHE' : 'CA HA CB HB1 HB2 CG CD1 HD1 CD2 HD2 CE1 HE1 CE2 HE2 CZ HZ',
    'PRO' : 'CA HA CB HB1 HB2 CD HD1 HD2 CG HG1 HG2',
    'SER' : 'CA HA CB HB1 HB2',
    'THR' : 'CA HA CB HB CG2 HG21 HG22 HG23',
    'TRP' : 'CA HA CB HB1 HB2 CG CD1 HD1 CD2 CE3 HE3 CZ3 HZ3 CH2 HH2 CZ2 HZ2',
    'TYR' : 'CA HA CB HB1 HB2 CG CD1 HD1 CD2 HD2 CE1 HE1',
    'VAL' : 'CA HA CB HB CG1 HG11 HG12 HG13 CG2 HG21 HG22 HG23',

    'POPC': 'C12 C11 H11A H11B C1 HA HB C2 HS C22 H2R H2S C3 HX HY C32 H2X H2Y C23 H3R H3S C24 H4R H4S C25 H5R H5S C26 H6R H6S C27 H7R H7S C28 H8R H8S C29 C210 C211 H11R H11S C212 H12R H12S C213 H13R H13S C214 H14R H14S C215 H15R H15S C216 H16R H16S C217 H17R H17S C218 H18R H18S H18T C33 H3X H3Y C34 H4X H4Y C35 H5X H5Y C36 H6X H6Y C37 H7X H7Y C38 H8X H8Y C39 H9X H9Y C310 H10X H10Y C311 H11X H11Y C312 H12X H12Y C313 H13X H13Y C314 H14X H14Y C315 H15X H15Y C316 H16X H16Y H16Z',
    'POPE': 'C12 H12A H12B C11 H11A H11B C1 HA HB C2 HS C22 H2R H2S C3 HX HY C32 H2X H2Y C23 H3R H3S C24 H4R H4S C25 H5R H5S C26 H6R H6S C27 H7R H7S C28 H8R H8S C29 C210 C211 H11R H11S C212 H12R H12S C213 H13R H13S C214 H14R H14S C215 H15R H15S C216 H16R H16S C217 H17R H17S C218 H18R H18S H18T C33 H3X H3Y C34 H4X H4Y C35 H5X H5Y C36 H6X H6Y C37 H7X H7Y C38 H8X H8Y C39 H9X H9Y C310 H10X H10Y C311 H11X H11Y C312 H12X H12Y C313 H13X H13Y C314 H14X H14Y C315 H15X H15Y C316 H16X H16Y H16Z',

    # DOPC & DOPE are the same
    'DOPC': 'C12 C11 H11A H11B C1 HA HB C2 HS C22 H2R H2S C3 HX HY C32 H2X H2Y C23 H3R H3S C24 H4R H4S C25 H5R H5S C26 H6R H6S C27 H7R H7S C28 H8R H8S C29 C210 C211 H11R H11S C212 H12R H12S C213 H13R H13S C214 H14R H14S C215 H15R H15S C216 H16R H16S C217 H17R H17S C218 H18R H18S H18T C33 H3X H3Y C34 H4X H4Y C35 H5X H5Y C36 H6X H6Y C37 H7X H7Y C38 H8X H8Y C39 C310 C311 H11X H11Y C312 H12X H12Y C313 H13X H13Y C314 H14X H14Y C315 H15X H15Y C316 H16X H16Y C317 H17X H17Y C318 H18X H18Y H18Z'  ,
    'DOPE': 'C12 C11 H11A H11B C1 HA HB C2 HS C22 H2R H2S C3 HX HY C32 H2X H2Y C23 H3R H3S C24 H4R H4S C25 H5R H5S C26 H6R H6S C27 H7R H7S C28 H8R H8S C29 C210 C211 H11R H11S C212 H12R H12S C213 H13R H13S C214 H14R H14S C215 H15R H15S C216 H16R H16S C217 H17R H17S C218 H18R H18S H18T C33 H3X H3Y C34 H4X H4Y C35 H5X H5Y C36 H6X H6Y C37 H7X H7Y C38 H8X H8Y C39 C310 C311 H11X H11Y C312 H12X H12Y C313 H13X H13Y C314 H14X H14Y C315 H15X H15Y C316 H16X H16Y C317 H17X H17Y C318 H18X H18Y H18Z',

    'POPI': 'C12 H2 C13 H3 C14 H4 C15 H5 C16 H6 C11 H1 C1 HA HB C2 HS C22 H2R H2S C3 HX HY C32 H2X H2Y C23 H3R H3S C24 H4R H4S C25 H5R H5S C26 H6R H6S C27 H7R H7S C28 H8R H8S C29 C210 C211 H11R H11S C212 H12R H12S C213 H13R H13S C214 H14R H14S C215 H15R H15S C216 H16R H16S C217 H17R H17S C218 H18R H18S H18T C33 H3X H3Y C34 H4X H4Y C35 H5X H5Y C36 H6X H6Y C37 H7X H7Y C38 H8X H8Y C39 H9X H9Y C310 H10X H10Y C311 H11X H11Y C312 H12X H12Y C313 H13X H13Y C314 H14X H14Y C315 H15X H15Y C316 H16X H16Y H16Z',
    
    'POPS': 'C12 C11 H11A H11B C1 HA HB C2 HS C22 H2R H2S C3 HX HY C32 H2X H2Y C23 H3R H3S C24 H4R H4S C25 H5R H5S C26 H6R H6S C27 H7R H7S C28 H8R H8S C29 C210 C211 H11R H11S C212 H12R H12S C213 H13R H13S C214 H14R H14S C215 H15R H15S C216 H16R H16S C217 H17R H17S C218 H18R H18S H18T C33 H3X H3Y C34 H4X H4Y C35 H5X H5Y C36 H6X H6Y C37 H7X H7Y C38 H8X H8Y C39 H9X H9Y C310 H10X H10Y C311 H11X H11Y C312 H12X H12Y C313 H13X H13Y C314 H14X H14Y C315 H15X H15Y C316 H16X H16Y H16Z',
    
    'CER180': 'H1S H1T H2S C4S C5S C6S H6S H6T C7S H7S H7T C8S H8S H8T C9S H9S H9T C10S H10S H10T C11S H11S H11T C12S H12S H12T C13S H13S H13T C14S H14S H14T C15S H15S H15T C16S H16S H16T C17S H17S H17T C18S H18S H18T H18U C2F H2F H2G C3F H3F H3G C4F H4F H4G C5F H5F H5G C6F H6F H6G C7F H7F H7G C8F H8F H8G C9F H9F H9G C10F H10F H10G C11F H11F H11G C12F H12F H12G C13F H13F H13G C14F H14F H14G C15F H15F H15G C16F H16F H16G C17F H17F H17G C18F H18F H18G H18H',
    
    'CHL1': 'C3 H3 C4 H4A H4B C5 C6 C7 H7A H7B C8 H8 C14 H14 C15 H15A H15B C16 H16A H16B C17 H17 C13 C18 H18A H18B H18C C12 H12A H12B C11 H11A H11B C9 H9 C10 C19 H19A H19B H19C C1 H1A H1B C2 H2A H2B C20 H20 C21 H21A H21B H21C C22 H22A H22B C23 H23A H23B C24 H24A H24B C25 H25 C26 H26A H26B H26C C27 H27A H27B H27C',
    
    # provided by Reza    
    'PSM': 'C11 H11A H11B C12 H12A H12B C13 H13A H13B H13C C14 H14A H14B H14C C15 H15A H15B H15C C1S H1S H1T C2S H2S C3 H3S C2F H2F H2G C3F H3F H3G C4F H4F H4G C5F H5F H5G C6F H6F H6G C7F H7F H7G C8F H8F H8G C9F H9F H9G C10F H10F H10G C11F H11F H11G C12F H12F H12G C13F H13F H13G C14F H14F H14G C15F H15F H15G C16F H16F H16G H16H C4S H4S C5S H5S C6S H6S H6T C7S H7S H7T C8S H8S H8T C9S H9S H9T C10S H10S H10T C11S H11S H11T C12S H12S H12T C13S H13S H13T C14S H14S H14T C15S H15S H15T C16S H16S H16T C17S H17S H17T C18S H18S H18T H18U',
    'SAPI24': 'C11 H1 C12 H2 C13 H3 C14 H4 C15 H5 C16 H6 C1 HA HB C2 H3 C3 HX HY C22 H2S H2R C23 H3S H3R C24 H4S H4R C25 H5R C26 H6R C27 H7S H7R C28 H8R C29 H9R C210 H10S H10R C211 H11R C212 H12R C213 H13S H13R C214 H14R C215 H15R C216 H16S H16R C217 H17S H17R C218 H18S H18R C219 H19S H19R C220 H20S H20R H20T C32 H2X H2Y C33 H3X H3Y C34 H4X H4Y C35 H5X H5Y C36 H6X H6Y C37 H7X H7Y C38 H8X H8Y C39 H9X H9Y C310 H10X H10Y C311 H11X H11Y C312 H12X H12Y C313 H13X H13Y C314 H14X H14Y C315 H15X H15Y C316 H16X H16Y C317 H17X H17Y C318 H18X H18Y H18Z',

    ## Not used/no head list provided
    'CHOL': 'C3 H3 C4 H4A H4B C5 C6 C7 H7A H7B C8 H8 C14 H14 C15 H15A H15B C16 H16A H16B C17 H17 C13 C18 H18A H18B H18C C12 H12A H12B C11 H11A H11B C9 H9 C10 C19 H19A H19B H19C C1 H1A H1B C2 H2A H2B C20 H20 C21 H21A H21B H21C C22 H22A H22B C23 H23A H23B C24 H24A H24B C25 H25 C26 H26A H26B H26C C27 H27A H27B H27C',
    
    ### The following head atoms will be removed from their corresponding list above
    # e.g. when analyzing POPC, the selected atoms will be those in 'POPC' but NOT in 'POPC_HEAD'

    # *PC heads are the same
    # 'POPC_HEAD': 'C11 C12 C13 C14 C15 H11A H11B', # from old charmm script
    'POPC_HEAD': 'C12 C11 H11A H11B C1 HA HB C2 HS',
    'DOPC_HEAD': 'C12 C11 H11A H11B C1 HA HB C2 HS',

    # *PE heads are the same
    'POPE_HEAD': 'C12 H12A H12B C11 H11A H11B C1 HA HB C2 HS C3 HX HY',
    'DOPE_HEAD': 'C12 H12A H12B C11 H11A H11B C1 HA HB C2 HS C3 HX HY',

    'POPI_HEAD': 'C12 H2 C13 H3 C14 H4 C15 H5 C16 H6 C11 H1 C1 HA HB C2 HS C3 HX HY',
    'POPS_HEAD': 'C12 C11 H11A H11B C1 HA HB C2 HS  C3 HX HY',
    'CER180_HEAD': 'H1S H1T H2S C4S',
    'CHL1_HEAD': '', # no head to be analyzed  
    
    'PSM_HEAD': 'C11 H11A H11B C12 H12A H12B C13 H13A H13B H13C C14 H14A H14B H14C C15 H15A H15B H15C C1S H1S H1T C2S H2S C3 H3S',
    'SAPI24_HEAD':'C11 H1 C12 H2 C13 H3 C14 H4 C15 H5 C16 H6 C1 HA HB C2 H3 C3 HX HY',

    # Not used
    # 'POPG_HEAD': 'C11 C12 C13 H11A H11B H12A H13A H13B', # from old charmm script
}
