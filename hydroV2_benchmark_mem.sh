#!/bin/bash

## For testing hydrophobic analysis with regard to different sizes of memory 

# Remeory to create local config file first
CONFIG="local/10582dcd.json"  # used phospho_d_pure_popc_read.dcd, 10582 frames, 13.9GB
PY_BIN="./hydrophobic_ana_v2.py"

REPORT=hydro_MEM_performance_report_$( hostname ).txt
echo "============ ${PY_BIN} ============" > $REPORT

for mem in 1 2 4 8  # test with different sizes of memory
do
    echo "Running with $mem GB ..."
    echo "------------- MEM: $mem GB  --------------" >> $REPORT
    sed "s/\"mem_limit\": 32,/\"mem_limit\": ${mem},/g" $CONFIG > $CONFIG'_sed'
    { time $PY_BIN $CONFIG'_sed' >> $REPORT ;} 2>> $REPORT

    echo -e "----------------------------------------------------\n" >> $REPORT
done

