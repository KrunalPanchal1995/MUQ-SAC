#!/bin/bash
# Define Chemkin running environment
#Uncomment the below line if not added to .bashrc file
#. /home/krithika/Desktop/ANSYS/ansys_inc/v231/reaction/chemkinpro.linuxx8664/bin/chemkinpro_setup.ksh

# Delete any intermediate files left over from a previous run
rm -f /home/krithika/MUQ-SAC/ISSUES_IN_CHEMKIN/CASE_6_rev/output/XMLdata_H2_0.zip

rm -f /home/krithika/MUQ-SAC/ISSUES_IN_CHEMKIN/CASE_6_rev/output/chemkindata.dtd

# Copy files to working directory
ln -s /home/krithika/Desktop/ANSYS/ansys_inc/v231/reaction/chemkinpro.linuxx8664/data/chemkindata.dtd /home/krithika/MUQ-SAC/ISSUES_IN_CHEMKIN/CASE_6_rev/output/chemkindata.dtd
# Run the job
CHEMKIN_MODE=Pro
export CHEMKIN_MODE
CKReactorGenericClosed -i /home/krithika/MUQ-SAC/ISSUES_IN_CHEMKIN/CASE_6_rev/output/H2.inp -o /home/krithika/MUQ-SAC/ISSUES_IN_CHEMKIN/CASE_6_rev/output/H2_run_0.out -x /home/krithika/MUQ-SAC/ISSUES_IN_CHEMKIN/CASE_6_rev/output/XMLdata_H2_0.zip Pro -m H2.out.mon -c /home/krithika/MUQ-SAC/ISSUES_IN_CHEMKIN/CASE_6_rev/output/H2_gas.asc
