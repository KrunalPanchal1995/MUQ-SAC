#!/bin/bash
# Define Chemkin running environment
#Uncomment the below line if not added to .bashrc file
#. /home/krithika/Desktop/ANSYS/ansys_inc/v231/reaction/chemkinpro.linuxx8664/bin/chemkinpro_setup.ksh

# Delete any intermediate files left over from a previous run
rm -f /data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/PRS_INVESTIGATION_ID_1029_TIG_CHEMKIN/Opt/case-6/0/output/XMLdata_H2_0.zip

rm -f /data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/PRS_INVESTIGATION_ID_1029_TIG_CHEMKIN/Opt/case-6/0/output/chemkindata.dtd

# Copy files to working directory
ln -s /home/krithika/Desktop/ANSYS/ansys_inc/v231/reaction/chemkinpro.linuxx8664/data/chemkindata.dtd /data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/PRS_INVESTIGATION_ID_1029_TIG_CHEMKIN/Opt/case-6/0/output/chemkindata.dtd
# Run the job
CHEMKIN_MODE=Pro
export CHEMKIN_MODE
CKReactorGenericClosed -i /data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/PRS_INVESTIGATION_ID_1029_TIG_CHEMKIN/Opt/case-6/0/output/H2.inp -o /data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/PRS_INVESTIGATION_ID_1029_TIG_CHEMKIN/Opt/case-6/0/output/H2_run_0.out -x /data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/PRS_INVESTIGATION_ID_1029_TIG_CHEMKIN/Opt/case-6/0/output/XMLdata_H2_0.zip Pro -m H2.out.mon -c /data/TEST-THERMO-sens/FFCM_THERMO_UQ_OPT/PRS_INVESTIGATION_ID_1029_TIG_CHEMKIN/Opt/case-6/0/output/H2_gas.asc
