#!/bin/bash
# Define Chemkin running environment
#Uncomment the below line if not added to .bashrc file
#. /home/krithika/Desktop/ANSYS/ansys_inc/v231/reaction/chemkinpro.linuxx8664/bin/chemkinpro_setup.ksh

# Extract solution data to CKCSV
GetSolution -nosen -norop -all -p CKSolnList.txt /home/krithika/MUQ-SAC/ISSUES_IN_CHEMKIN/CASE_6_rev/output/XMLdata_H2_0.zip

# Delete postprocessor log file
rm -f /home/krithika/MUQ-SAC/ISSUES_IN_CHEMKIN/CASE_6_rev/output/ckpp_*.log

