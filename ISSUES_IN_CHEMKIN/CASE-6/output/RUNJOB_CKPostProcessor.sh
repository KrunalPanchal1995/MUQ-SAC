#!/bin/bash
# Define Chemkin running environment
#Uncomment the below line if not added to .bashrc file
#. /home/krithika/Desktop/ANSYS/ansys_inc/v231/reaction/chemkinpro.linuxx8664/bin/chemkinpro_setup.ksh

# Extract solution data to CKCSV
GetSolution -nosen -norop -all -p CKSolnList.txt /data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/PRS_INVESTIGATION_ID_1029_TIG_CHEMKIN/Opt/case-6/0/output/XMLdata_H2_0.zip

# Delete postprocessor log file
rm -f /data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/PRS_INVESTIGATION_ID_1029_TIG_CHEMKIN/Opt/case-6/0/output/ckpp_*.log

