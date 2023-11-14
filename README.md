# User Guide
## Basic usage:
Step 1. [Set up conda environment](./docs/notes.md#1-python-environment-setup)
Step 2. Prepare a config file
    - the config file sets parameters for the analysis; check [detailed explanation of each parameter](./config_template.json), 
    - [example config file](./config_example.json)
Step 3. Run one of the analysis scripts through command line: `hbond_ana.py`, `hydrophobic_ana.py`, or `cation_pi_ana.py`
```bash
(md-ana-py37) $ ./<analysis_code>  <configuration_file>

## for example, run hbond analysis
(md-ana-py37) $ ./hbond_ana.py config_example.json
```

## Addtional notes:
1. [Conda environment setup](./docs/notes.md#1-python-environment-setup)
2. [Performance tests](./docs/notes.md#2-performances)