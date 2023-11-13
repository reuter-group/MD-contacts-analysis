#!/bin/bash

PSF="test_data/phosphod_bound_to_pure_popc.old.psf"

f100="test_data/100frames.dcd"
f500="test_data/500frames.dcd"
### On my laptop
# f2501="local/last_50ns.dcd"
# f10582="local/phospho_d_pure_popc_read.dcd"

### On Betzy
f2501="/cluster/projects/nn4700k/dandan/md/last_50ns.dcd"
f10582="/cluster/projects/nn4700k/dandan/md/phospho_d_pure_popc_read.dcd"


for PY_BIN in "./cation_pi_ana.py" "./hydrophobic_ana.py"
do 
    REPORT=performance_report_$( hostname )_$( basename ${PY_BIN} ).txt
    echo "============ ${PY_BIN} ============" > $REPORT

    if [ "${PY_BIN}" == "./cation_pi_ana.py" ]; then
        RES="test_data/cation_pi_candidates.txt"
    else
        RES="test_data/hydroph_candidates.txt"
    fi
    echo "Running $PY_BIN ..."
    for input in $f100 $f500 $f2501 $f10582
    do
        echo "Running on $input ..."
        echo "------------- $input --------------" >> $REPORT
        { time $PY_BIN $PSF $input $RES 2 >> $REPORT ;} 2>> $REPORT

        echo -e "----------------------------------------------------\n" >> $REPORT
    done
done