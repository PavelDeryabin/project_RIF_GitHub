#!/bin/bash 

# GEDIT deconvolution algorithm was downloaded from GitHub repository https://github.com/BNadel/GEDITExpanded
# use python3 and run GEDITExpanded-main/GEDITv3.0/GEDIT3.py script as follows

python3 GEDIT3.py \
-mix [PATH_TO_MATRIX_TO_BE_DECONVOLUTED] \
-ref [PATH_TO_processed.data_DIRECTORY]/seu.PRJNA1141235.fig3.reference.sliced.csv \
-outFile [PATH_TO_processed.data_DIRECTORY]/dcnv.results
