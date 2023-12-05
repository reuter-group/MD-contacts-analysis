# User Guide

### Prerequisites
- For installing the scripts, [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) or Anaconda is required. 

### Step 1. Set up running environment
1.1 Download the source code and create the conda environment 
```
git clone https://github.com/reuter-group/MD-contacts-analysis.git
cd /path/to/MD-contacts-analysis
conda create -f md-ana-py37.yml
```
1.2 Activate the environment
```
conda activate md-ana-py37
```

### Step 2. Prepare a config file
The config file sets parameters for the analysis; check [detailed explanation of each parameter](./config_template.json). An example config file `config_example.json`:
```
{
    "PSF_FILE": "test_data/reza_data.psf",

    "DCD_FILE": [
        "test_data/reza_data_64frames.dcd"
    ],

    "RESIDUE_FILE": "test_data/reza_data_short_candidates.txt",

    "hbond_lipids": ["POPC", "POPE", "POPI", "POPS", "CHL1", "PSM", "SAPI24"],
    "hydrophobic_lipids": ["POPC", "POPE", "CHL1", "PSM", "SAPI24"],
    "cation_pi_lipids": ["POPC", "POPE", "POPS", "SAPI24", "PSM"],    

    "lifetime": 2, 
    
    "mem_limit": 4,

    "start": 0,
    "stop": null, 
    "step": 1,

    "output_file_prefix": "test-reza"
}
```

### Step 3. Run analysis script 
Choose from one of the following scripts to run:
- Hydrogen bonds: `hbond_ana.py`
- Hydrophobic contacts: `hydrophobic_ana.py`
- Cation-Pi interactions: `cation_pi_ana.py`

```bash
(md-ana-py37) $ ./<analysis_code>  <configuration_file>

## for example, run hbond analysis
(md-ana-py37) $ ./hbond_ana.py config_example.json
```
