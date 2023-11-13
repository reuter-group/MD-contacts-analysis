# User Guide
## Basic usage:
1. [Set up the environment](https://git.app.uib.no/reuter-group/md_analysis/-/blob/master/docs/notes.md#1-python-environment-setup)
2. Prepare a configuration file
    - A _configuration file_ is for setting parameters for the analysis; check [detailed explanation of each parameter](https://git.app.uib.no/reuter-group/md_analysis/-/blob/master/config_template.json), 
    - NOTE: the configuration file has a **strict** format JSON, so don't modify anything that you are not sure about
    - when you get familar with the parameters, you can use a shorter example [like this](https://git.app.uib.no/reuter-group/md_analysis/-/blob/master/config_example.json)
3. Run analysis code in command line, where the _analysis code_ is one of the three analysis scripts: `hbond_ana.py`, `hydrophobic_ana.py`, or `cation_pi_ana.py`
```bash
(md-ana-py37) $ ./<analysis_code>  <configuration_file>

## for example, run hbond analysis
(md-ana-py37) $ ./hbond_ana.py config_example.json
```


## Addtional notes:
1. [Environment setup](https://git.app.uib.no/reuter-group/md_analysis/-/blob/master/docs/notes.md#1-python-environment-setup)
2. [Performance tests](https://git.app.uib.no/reuter-group/md_analysis/-/blob/master/docs/notes.md#2-performances)