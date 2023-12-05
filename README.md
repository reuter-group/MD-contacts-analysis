# User Guide

### Step 1. Set up running environment

- Download the source code and create the environment with `conda`
```
conda create -f md-ana-py37.yml

```
- Activate the environment
```
$ conda activate md-ana-py37
```


### Step 2. Prepare a config file
- the config file sets parameters for the analysis; check [detailed explanation of each parameter](./config_template.json)
- an example config file `config_example.json`
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

### Step 3. Run one of the analysis scripts in command line: `hbond_ana.py`, `hydrophobic_ana.py`, or `cation_pi_ana.py`
```bash
(md-ana-py37) $ ./<analysis_code>  <configuration_file>

## for example, run hbond analysis
(md-ana-py37) $ ./hbond_ana.py config_example.json
```
