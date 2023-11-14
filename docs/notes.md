
## Addtional notes:
1. [Environment setup](#1-python-environment-setup)
2. [Performance tests](#2-performances)
    

### 1. Python environment setup
- create the environment with `conda`
```
$ conda env create -f md-ana-py37.yml  # on Mac

# or on Linux
conda create -f md-ana-py37.yml

```
- activate the environment
```
$ conda activate md-ana-py37
```


### 2. Performances 

#### 2.1 Inputs

- DCD files:
    - [100 frames](./test_data/100frames.dcd): 132 MB
    - [500 frames](./test_data/500frames.dcd): 657 MB
    - 2501 frames: 3.3 GB
    - 10582 frames: 13.9 GB

Every test was running the analysis with selecting 2 consecutive frames. The running time was reported by the Unix `time` command.
Example:

```
$ time ./cation_pi_ana.py test_data/phosphod_bound_to_pure_popc.old.psf test_data/500frames.dcd test_data/cation_pi_candidates.txt 2
```
For more detials, check [test script `benchmark.sh`](./benchmark.sh)


#### 2.2 Platforms
* For the [old scripts](./charmm_scripts): `charmm`+`csh`(`sed`, `awk`), tested on `doggpil.cbu.uib.no` where the required charmm is installed
* The new version (this repository) was tested on two platforms:
    * Laptop: MacBook Pro, 4 cores (Intel Core i5), 8 GB memory
    * HPC: one `Betzy` compute node, 128 cores, 256 GB memory, used `--cpus-per-task=128`



#### 2.3 Hydrophobic Contacts Analysis

| Input       |  Old script   |  Laptop   |  HPC  | 
|-------------|:-------------:|:---------:|:-----:|
| 100 frames  | 87            | 37        |  35  |
| 500 frames  | 238           | 40        |  41  |
| 2501 frames | 632           | 94         |  46  |
| 10582 frames| 2846          | (out of memory)* | 68  |

- time unit: second
- [More details about HPC's result](./performance/performance_report_betzy.b1184_hydrophobic.txt)
- Note: 
    - the results of these two versions are slightly different ( < 1% ), mainly due to different definitions of lifetime cutoff 
    - [hydrophobic candidates used for this test](https://git.app.uib.no/reuter-group/md_analysis/-/blob/4362a1e6074ed818c36b692129fc94fd83b423db/hydrop_candi.py)
- *in newer version (since 0596e69244baeaa213460615a8fe5f9683106ddf) where specifying memory limit is supported, the performance is 463s (7~8m). 


#### 2.4 Cation-pi Contacts Analysis

| Input       |  Old script   |  Laptop   |  HPC  |
|-------------|:-------------:|:---------:|:-----:|
| 100 frames  |  619          | 10        |  10  |
| 500 frames  |  1098         | 11        |  10  |
| 2501 frames |  3692         | 12        |  11  |
| 10582 frames|  ? ( >10000 )  | 21        |  17  |

- time unit: second
- [More details about HPC's result](./performance/performance_report_betzy.b1184_cation_pi.txt)
- Note: the results of two versions are different ( <10%, for 2501 frames' case), for two reasons:
    - the old script pre-selected lipid candidates from the first 1% number of frames, which leads to missing candidates in most cases.
    - the old script may not have the lifetime cutoff processing


#### 2.5 Hydrogen Bond Analysis

| Input       |  Old script   |  Laptop   | 
|-------------|:-------------:|:---------:|
| 100 frames  |  393          | 6        |
| 500 frames  |  424         | 8        |
| 2501 frames |  581         | 21        |
| 10582 frames|  1202        | 100       |

- time unit: second
- [More details about Laptop's result](./performance/performance_report_Laptop_hbond_ana.py.txt)
- the code implementation is based on [MDAnalysis's HBond module](https://www.mdanalysis.org/docs/documentation_pages/analysis/hbond_analysis.html)

