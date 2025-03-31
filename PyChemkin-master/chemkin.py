#!/usr/bin/env python
# encoding: utf-8

"""
This module contains several helper functions useful when running Chemkin jobs
from the command line.
"""
import os
import stat
import os.path
import subprocess
import csv
import numpy
import logging
import sys
import numpy as np

# The directory in which Chemkin is installed
chemkin_path = os.path.abspath(os.path.join(os.path.dirname(__file__),'chemkin_path'))

if os.path.exists(chemkin_path):
    with open(chemkin_path, 'r') as f:
        CHEMKIN_DIR = os.path.abspath(f.readline()).strip("\n")


os.environ["CHEMKIN_BIN"] = os.path.join(CHEMKIN_DIR,"bin")
os.environ["LD_LIBRARY_PATH"] = "/home/krithika/Desktop/ANSYS/ansys_inc/v231/tp/IntelMKL/2021.3.0/linx64/lib/intel64"
os.environ["LD_LIBRARY_PATH"] = "/home/krithika/Desktop/ANSYS/ansys_inc/v231/tp/IntelCompiler/2019.3.199/linx64/lib/intel64"
os.environ["LD_LIBRARY_PATH"] = "/home/krithika/Desktop/ANSYS/LIB_FILES_FOR_CHEMKIN_API"
# Path to MKL libraries, update as needed
chemkin_pro_setup = os.path.join(CHEMKIN_DIR, "bin", "chemkinpro_setup.ksh")

#print("Step 1: Setting up the chemkin environment\n")
# Run the setup script to configure Chemkin Pro
subprocess.call(f"bash -c 'source {chemkin_pro_setup}'", shell=True)
#print("Step 1: Executed\n")
################################################################################


# The preamble to each Chemkin execution shell script
CHEMKIN_SCRIPT_PREAMBLE = """#!/bin/bash
# Define Chemkin running environment
#Uncomment the below line if not added to .bashrc file
#. {0}

""".format(os.path.join(CHEMKIN_DIR, 'bin', 'chemkinpro_setup.ksh'))

################################################################################

# Set logging parameters

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Create console handler and set level to debug; send everything to stdout
# rather than stderr
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)

logging.addLevelName(logging.CRITICAL, 'Critical: ')
logging.addLevelName(logging.ERROR, 'Error: ')
logging.addLevelName(logging.WARNING, 'Warning: ')
logging.addLevelName(logging.INFO, '')

formatter = logging.Formatter('%(levelname)s%(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

################################################################################

class ChemkinJob(object):
    """
    A single (gas-phase) Chemkin job.
    """
    
    def __init__(self, name, chemFile, tempDir):
        self.name = name
        self.case =  0
        self.chemFile = os.path.abspath(chemFile)
        self.tempDir = os.path.abspath(tempDir).strip("\n")
        if not os.path.exists(self.tempDir):
            os.makedirs(self.tempDir)
        
    @property
    def ascFile(self):
        return os.path.join(self.tempDir, '{0}_gas.asc'.format(self.name))
    
    @property
    def inputFile(self):
        return os.path.join(self.tempDir, '{0}.inp'.format(self.name))
    
    @property
    def outputFile(self):
        return os.path.join(self.tempDir, '{0}_run_{1}.out'.format(self.name,self.case))
    
    @property
    def monFile(self):
        return os.path.join('{0}.out.mon'.format(self.name))
    
    @property
    def ckcsvFile(self):
        return os.path.join(self.tempDir, 'CKSoln.ckcsv'.format(self.case))

    @property
    def preferenceFile(self):
        return os.path.join(self.tempDir, 'CKSolnList.txt')
    
    @property
    def dataZipFile(self):
        return os.path.join(self.tempDir, 'XMLdata_{0}_{1}.zip'.format(self.name,self.case))    
    
    def preprocess(self):
        """
        Run the Chemkin preprocessor on the chemistry input file.
        """
        # Write the preprocessor execution script to a file
        scriptFile = os.path.join(self.tempDir, 'RUNJOB_CKPreProcess.sh')
        with open(scriptFile, 'w') as stream:
            stream.write(CHEMKIN_SCRIPT_PREAMBLE)
            stream.write('{0} -i {1} -o {2}_gas.out -c {2}_gas.asc\n\n'.format(
                os.path.join(CHEMKIN_DIR, 'bin', 'chem'),
                os.path.abspath(self.chemFile),
                os.path.join(self.tempDir, self.name),
            ))
        os.chmod(scriptFile, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
        # Execute the preprocessor script
        #print("RUn file created")
        process = subprocess.Popen(('/bin/bash', scriptFile), cwd=self.tempDir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print(stdout)
            print(stderr)
            quit()
        
    def run(self, jobParams, model, pro=True):
        """
        Run a Chemkin job.
        """
        dtdFileSrc = os.path.join(CHEMKIN_DIR, 'data', 'chemkindata.dtd')
        dtdFileDest = os.path.join(self.tempDir, 'chemkindata.dtd')
        scriptFile = os.path.join(self.tempDir, 'RUNJOB_CKRunProcessor.sh')
        # Write the job parameters to a file
        with open(self.inputFile, 'w') as stream:
            stream.write(jobParams)
        # Write the job execution script to a file
        with open(scriptFile, 'w') as stream:
            stream.write(CHEMKIN_SCRIPT_PREAMBLE)
            stream.write('# Delete any intermediate files left over from a previous run\n')
            stream.write('rm -f {0}\n\n'.format(self.dataZipFile))
            stream.write('rm -f {0}\n\n'.format(dtdFileDest))
            stream.write('# Copy files to working directory\n')
            stream.write('ln -s {0} {1}\n'.format(dtdFileSrc, dtdFileDest))
            stream.write('# Run the job\n')
            if pro:
                stream.write('CHEMKIN_MODE=Pro\n')
                stream.write('export CHEMKIN_MODE\n')
                stream.write('{0} -i {1} -o {2} -x {3} Pro -m {4} -c {5}\n'.format(model, self.inputFile, self.outputFile, self.dataZipFile, self.monFile, self.ascFile))
            else:
                stream.write('{0} -i {1} -o {2} -x {3} -m {4} -c {5}\n'.format(model, self.inputFile, self.outputFile, self.dataZipFile, self.monFile, self.ascFile))
        # Execute the job script
        subprocess.call(('/bin/bash', scriptFile), cwd=self.tempDir)
    def updateCase(self):
        self.case += 1

    def postprocess(self, sens=False, rop=False, all=False, transpose=True):
        """
        Run the Chemkin postprocessor.
        """
        # Set default preferences and write them to the preference file
        with open(self.preferenceFile, 'w') as stream:
            stream.write('VARIABLE VAR ALL\n')
            stream.write('VARIABLE SEN ALL\n')
            stream.write('VARIABLE ROP ALL\n') 
            stream.write('VARIABLE volume  1  0  0\n')
            stream.write('VARIABLE external_surface_area  0  0  0\n')
            stream.write('VARIABLE volumetric_heat_production_rate  0  0  0\n')
            stream.write('VARIABLE surface_temperature  0  0  0\n')
            stream.write('VARIABLE heat_loss_rate  0  0  0\n')
            stream.write('VARIABLE temperature  1  0  0\n')
            stream.write('VARIABLE mass  0  0  0\n')
            stream.write('VARIABLE pressure  1  0  0\n')
            stream.write('VARIABLE molar_conversion  0  0  0\n')
            stream.write('VARIABLE net_reaction_rate  0  0  0\n')
            stream.write('VARIABLE heat_production_per_reaction  0  0  0\n')
            stream.write('VARIABLE molecular_weight  0  0  0\n')
            stream.write('VARIABLE mass_density  0  0  0\n')
            stream.write('VARIABLE mixture_enthalpy  0  0  0\n')
            stream.write('VARIABLE sensible_enthalpy  0  0  0\n')
            stream.write('VARIABLE unburned_hydrocarbons  0  0  0\n')
            stream.write('VARIABLE volatile_organic_compounds  0  0  0\n')
            stream.write('VARIABLE exhaust_in_ppmvd  0  0  0\n')
            stream.write('VARIABLE all_single_point_values  0  0  0\n')
            stream.write('VARIABLE use_reaction_string_instead_of_index  1  0  0\n')
            stream.write('UNIT  Time  (sec)\n')
            stream.write('UNIT  Temperature  (K)\n')
            stream.write('UNIT  Pressure  (bar)\n')
            stream.write('UNIT  Volume  (cm3)\n')

        # Write the postprocessor execution script to a file
        scriptFile = os.path.join(self.tempDir, 'RUNJOB_CKPostProcessor.sh')
        with open(scriptFile, 'w') as stream:
            stream.write(CHEMKIN_SCRIPT_PREAMBLE)
            stream.write('# Extract solution data to CKCSV\n')
            stream.write('GetSolution {1}{2}{3}{4}{0}\n\n'.format(
                self.dataZipFile, 
                '-nosen ' if not sens else '',
                '-norop ' if not rop else '',
                '-all ' if all else '',
                '-p'+' CKSolnList.txt ',
            ))
            if transpose:
                stream.write('# Transpose the data to CSV\n')
                stream.write('CKSolnTranspose -column 500 {0}\n\n'.format(self.ckcsvFile))
            stream.write('# Delete postprocessor log file\n')
            stream.write('rm -f {0}\n\n'.format(os.path.join(self.tempDir, 'ckpp_*.log')))
        
        # Execute the postprocessor script
        subprocess.call(('/bin/bash', scriptFile), cwd=self.tempDir)
        
    def writeInputHomogeneousBatch(self,problemType, reactants, temperature, pressure, endTime, 
                      Continuations=False, typeContinuation = None, Tlist = [], Plist = [],
                      variableVolume=False, variableVolumeProfile = None, 
                      solverTimeStepProfile = None,variableVolumeProfileType=None):
        """
        Write input file for homogeneous batch reactor
        """
        
        # Problem definition block 
                    
        input_stream=("""! 
! problem type definition
!""")

        if problemType.lower() == 'constrainVandSolveE'.lower():
            input_stream+=("""
CONV   ! Constrain Volume And Solve Energy Equation
ENRG   ! Solve Gas Energy Equation""")
        elif problemType.lower() == 'constrainPandSolveE'.lower():
            input_stream+=("""
ENRG   ! Solve Gas Energy Equation""")
        elif problemType.lower() == 'constrainPandT'.lower():
            input_stream+=("""
TGIV   ! Constrain Pressure And Temperature
""")
        elif problemType.lower() == 'constrainVandT'.lower():
            input_stream+=("""
COTV   ! Constrain Volume And Temperature
TGIV   ! Fix Gas Temperature""")
             
        # Solver type definition block 
        
        input_stream+=("""
TRAN   ! Transient Solver""")
        
        # Physical property
        
        input_stream+=("""
! 
! physical property
! 
!Surface_Temperature   ! Surface Temperature Same as Gas Temperature
IFAC 0.1   ! Ignition Noise Filtering Factor
""")
            
        input_stream+=('PRES {0:g}   ! Pressure (atm)\n'.format(pressure/101325)) #Make sure the pressure units are matching with the exp conditions
        input_stream+=('QLOS 0.0   ! Heat Loss (cal/sec)\n')
        input_stream+=('TEMP {0:.2f}   ! Temperature (K)\n'.format(temperature))
        input_stream+=("""VOL 1.0   ! Volume (cm3)""")
        if variableVolume:
            if variableVolumeProfileType == None:
                with open(variableVolumeProfile, 'rb') as csvfile:
                    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                    for row in reader:
                        time = float(row[0].split(',')[0]) # (sec)
                        vol = float(row[0].split(',')[1]) # (cm3)
                        input_stream+=("""
VPRO {0:g} {1:g}   ! Volume (cm3)""".format(time,vol))
            elif variableVolumeProfileType == "ck_file":
                 with open(variableVolumeProfile, 'r') as csvfile:
                    reader = csvfile.readlines()
                    for row in reader:
                        #time = float(row[0].split(',')[0]) # (sec)
                        #vol = float(row[0].split(',')[1]) # (cm3)
                        input_stream+=("""
{0}   ! Volume (cm3)""".format(row))
            else:
                time = variableVolumeProfile.time
                volume = variableVolumeProfile.volume
                for i,t in enumerate(time):
                    input_stream+=("""
VPRO {0:g} {1:g}   ! Volume (cm3)""".format(t,volume[i]))
        
                       
        
        # Species property block
        
        input_stream+=("""
!
! species property
!
""")
        
        for reac , conc in reactants:            
            input_stream+=('REAC {0} {1:g} ! Reactant Fraction (mole fraction) \n'.format(reac,conc))
            
        # Solver control block
        
        input_stream+=("""! 
! solver control
! 
ADAP   ! Save Additional Adaptive Points
NNEG   !Forced Non-Negative Solutions
!ASTEPS 20   ! Use Solver Integration Steps (Uncomment this only when needed)
!ATLS 1.0E-6   ! Sensitivity Absolute Tolerance
!ATOL 1.0E-8   ! Absolute Tolerance
!MAXIT 4   ! Maximum Number of Iterations
!RTLS 0.0001   ! Sensitivity Relative Tolerance
!RTOL 1.0E-4   ! Relative Tolerance
!STPT 1.0E-6   ! Solver Maximum Step Time (sec)""")
        
        if solverTimeStepProfile:
            with open(solverTimeStepProfile, 'rb') as csvfile2:
                timeStepReader = csv.reader(csvfile2, delimiter=' ', quotechar='|')
                for row in timeStepReader:
                    time = float(row[0].split(',')[0]) # (sec)
                    vol = float(row[0].split(',')[1]) # (sec)
                    input_stream+=("""
STPTPRO {0:g} {1:g}           ! Solver Maximum Step Time (sec)""".format(time,vol))
            
        
        input_stream+=("""
!STPT {1:.7f}    ! Solver Maximum Step Time (sec)
TIME {0:.5f}                 ! End Time (sec)""".format(endTime,endTime/10000))
        
        input_stream+=("""
! 
! output control and other misc. property
! 
!CSVOUTPUT_KOUT CO2 FIXED_mass_fraction   ! Table Output Parameter
!CSVOUTPUT_KOUT H2O FIXED_mass_fraction   ! Table Output Parameter
!CSVOUTPUT_KOUT MB-C5H10O2 FIXED_mass_fraction   ! Table Output Parameter
!CSVOUTPUT_KOUT OH FIXED_mass_fraction   ! Table Output Parameter
!CSVOUTPUT_WRITE   ! Write CSV Output File
EPSR 0.01   ! Threshold for Rate of Production
EPSS 0.001   ! Threshold for Species Sensitivity
EPST 0.001   ! Threshold for Temperature Sensitivity
GFAC 1.0   ! Gas Reaction Rate Multiplier
IEEE_CHECK_USUAL 0   ! Disable Over/Underflow Check
PRNT 1   ! Print Level Control
!SDOUTPUT_WRITE   ! Write Standard Diagnostic Output
SIZE 100000000   ! Solution Data Block Size (bytes)
!XMLOUTPUT_WRITE   ! Write XML Solution File""")
              
        if Continuations:
            if numpy.array(Tlist).size:                
                for i in range(numpy.array(Tlist).shape[0]):
                    input_stream+=("""
{0}
END
TEMP {1:g}""".format(typeContinuation,numpy.array(Tlist)[i]))
            
            if numpy.array(Plist).size:         
                for i in range(numpy.array(Plist).shape[0]):
                    input_stream+=("""
{0}
END
PRES {1:g}""".format(typeContinuation,numpy.array(Plist)[i]/1.01325))
                                            
        input_stream+=('\nEND')
        
        return input_stream      

    def writeInputPlugFlow(self,problemType, reactants, 
                           startingAxialPosition, endingAxialPosition, diameter, momentum=True, massflowrate = None, sccmflowrate = None, 
                           temperature  = None, pressure  = None,  
                           temperatureProfile = None, pressureProfile = None):
    
        """
        Write input file for typical plug flow
        """
        
        # Problem definition block 
                    
        input_stream=("""
! 
! problem type definition
!
""")

        if momentum:
            input_stream+=("""
MOMEN ON   ! Turn on Momentum Equation
""")
        else:
            input_stream+=("""
MOMEN OFF   ! Turn off Momentum Equation
""")

        input_stream+=("""
PLUG   ! Plug Flow Reactor
RTIME ON   ! Turn on Residence Time Calculation
""")

        if problemType.lower() == 'FixGasTemperature'.lower():
            input_stream+=("""
TGIV   ! Fix Gas Temperature""")
        elif problemType.lower() == 'solveGasEnergy'.lower():
            input_stream+=("""
ENRG   ! Solve Gas Energy Equation""")
    
        
        # Physical property
        
        input_stream+=("""
! 
! physical property
! 
!Surface_Temperature   ! Surface Temperature Same as Gas Temperature
""")
        if massflowrate:
            input_stream+=('FLRT {0:g}   ! Mass Flow Rate (g/sec)\n'.format(massflowrate))
        elif sccmflowrate:
            input_stream+=('SCCM {0:g}   ! Volumetric Flow Rate in SCCM (standard-cm3/min@298.15K)\n'.format(sccmflowrate))
        else:
            raise Exception('Must indicate either a mass or sccm flow rate')
        
        if pressureProfile:
            with open(pressureProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    pos = float(row[0].split(',')[0]) # (cm)
                    pres = float(row[0].split(',')[1])*1.0/1.01325 # (atm)
                    input_stream+=("""
PPRO {0:g} {1:g}   ! Pressure (atm)""".format(pos,pres))
        else:
            input_stream+=('PRES {0:g}   ! Pressure (atm)\n'.format(pressure/1.01325))
            
        if temperatureProfile:
            with open(temperatureProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    pos = float(row[0].split(',')[0]) # (cm)
                    temp = float(row[0].split(',')[1]) # (K)
                    input_stream+=("""
TPRO {0:g} {1:g}   ! Temperature (K)""".format(pos,temp))
        else:
            input_stream+=('TEMP {0:g}   ! Temperature (K)'.format(temperature)) 
        
        # reactor dimension definition
        
        input_stream+=("""
VIS 0.0   ! Mixture Viscosity at Inlet (g/cm-sec)
! 
! reactor dimension definition
! 
""")
        input_stream+=('DIAM {0:g}   ! Diameter (cm)\n'.format(diameter))
        input_stream+=('XEND {0:g}   ! Ending Axial Position (cm)\n'.format(endingAxialPosition))
        input_stream+=('XSTR {0:g}   ! Starting Axial Position (cm)\n'.format(startingAxialPosition))
                       
        
        # Species property block
        
        input_stream+=("""
!
! species property
!
""")
        
        for reac , conc in reactants:            
            input_stream+=('REAC {0} {1:g} ! Reactant Fraction (mole fraction) \n'.format(reac,conc))
            
        # Solver control block
        
        input_stream+=("""! 
! 
! solver control
! 
ACHG 0.0   ! Maximum Absolute Change in Site Fractions
ADAP   ! Save Additional Adaptive Points
ATLS 1.0E-6   ! Sensitivity Absolute Tolerance
ATOL 1.0E-9   ! Absolute Tolerance
MAXIT 10   ! Maximum Number of Iterations
NNEG   ! Force Non-negative Solution
PSV 1.0E-8   ! Scaling Factor for Relaxing Surface Equations (cm/sec)
RCHG 1.0E-6   ! Maximum Relative Change in Site Fractions
RTLS 0.0001   ! Sensitivity Relative Tolerance
RTOL 1.0E-10   ! Relative Tolerance
TSTP 1.0   ! Initial Integration Step (cm""")          
        
        input_stream+=("""
! 
! output control and other misc. property
! 
GFAC 1.0   ! Gas Reaction Rate Multiplier""")
                                            
        input_stream+=('\nEND')
        
        return input_stream      
    
    def writeinputJSR(self,problemType, reactants, tau,endtime, volume, 
                           temperature  = None, pressure  = None,
                           Continuations=False, typeContinuation = None, Tlist = [], Plist = [],                           
                           temperatureProfile = None, pressureProfile = None):
    
        """
        Write input file for JSR
        """
        
        # Problem definition block 
                    
        input_stream=("""
! 
! problem type definition
!
TRAN   ! Transient Solver 
""")

        if problemType.lower() == 'FixGasTemperature'.lower():
            input_stream+=("""
TGIV   ! Fix Gas Temperature""")
        elif problemType.lower() == 'solveGasEnergy'.lower():
            input_stream+=("""
ENRG   ! Solve Gas Energy Equation""")
    
        
        # Physical property
        
        input_stream+=("""
! 
! physical property
! 
!Surface_Temperature   ! Surface Temperature Same as Gas Temperature
""")
                
        if pressureProfile:
            with open(pressureProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    pos = float(row[0].split(',')[0]) # (cm)
                    pres = float(row[0].split(',')[1])*1.0/1.01325 # (atm)
                    input_stream+=("""
PPRO {0:g} {1:g}   ! Pressure (atm)""".format(pos, pres))
        else:
            input_stream+=('PRES {0:g}   ! Pressure (atm)\n'.format(pressure/1.01325))
        
        # Residence time 
        input_stream+=('TAU {0:g}   ! Residence time (sec)\n'.format(tau))

        if temperatureProfile:
            with open(temperatureProfile, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for row in reader:
                    pos = float(row[0].split(',')[0]) # (cm)
                    temp = float(row[0].split(',')[1]) # (K)
                    input_stream+=("""
TPRO {0:g} {1:g}   ! Temperature (K)""".format(pos,temp))
        else:
            input_stream+=('TEMP {0:g}   ! Temperature (K)\n'.format(temperature)) 
        
        # volume
        input_stream+=('VOL {0:g}   ! Volume (cm3)\n'.format(volume))
        
        # Species property block
        
        input_stream+=("""
!
! species property
!
""")
        
        for reac , conc in reactants:            
            input_stream+=('REAC {0} {1:g} ! Reactant Fraction (mole fraction) \n'.format(reac,conc))
        # For transient solver you also need estimates of species initial gas fraction
        for reac , conc in reactants:            
            input_stream+=('XEST {0} {1:g} ! Initial Gas Fraction (mole fraction) \n'.format(reac,conc))
        
        # solver control    
        input_stream+=("""
! 
! solver control
! 
ADAP   ! Save Additional Adaptive Points
""")    
        input_stream+=('TIME {0:g}   ! End Time (sec) \n'.format(endtime))    
    
        # output control and other misc. property
        input_stream+=("""
! 
! output control and other misc. property
! 
GFAC 1.0   ! Gas Reaction Rate Multiplier
""")
                                            
        if Continuations:
            if numpy.array(Tlist).size:                
                for i in range(numpy.array(Tlist).shape[0]):
                    input_stream+=("""
{0}
END
TEMP {1:g}""".format(typeContinuation,numpy.array(Tlist)[i]))
            
            if numpy.array(Plist).size:         
                for i in range(numpy.array(Plist).shape[0]):
                    input_stream+=("""
{0}
END
PRES {1:g}""".format(typeContinuation,numpy.array(Plist)[i]/1.01325))
                                            
        input_stream+=('\nEND')        
        return input_stream      


################################################################################

def getVariable(ckcsvFile, variable=None):
    """
    Return the time and variable data (i.e. Pressure, Volume, or Temperature) from the given CKCSV file.  Returns the data
    in the form of [time_soln1, timesoln2, ...] [variable_data_soln1, variable_data_soln2,...]
    Returns time in units of [seconds]. Returns temperature in units of [K].
    Returns pressure in units of [bar]. Returns volume in units of [cm3].
    """
    if variable.lower() not in ['pressure','volume','temperature']:
        raise Exception('Can only parse Pressure, Volume, or Temperature from CKCSV file.')

    timeData = []; varData = []

    units = {'k': 1.0, 'bar': 1, 'mpa': 10, 'pa': 1e-5, 'atm': 1.01325, 'cm3': 1.0, 'm3': 1e6}
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            tokens = label.split('_')
            if tokens[0] == 'Time':
                tdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
                timeData.append(tdata)      

            if tokens[0].lower() == variable.lower():
                vdata = numpy.array([float(r) for r in row[2:]], numpy.float)
                vunits = row[1].strip()[1:-1].lower()
                vdata = vdata*units[vunits]
                varData.append(vdata)
    
    return timeData, varData

def getMoleFraction(ckcsvFile, species=[]):
    """
    Return the time and mole fraction data from the given CKCSV file.  
    Returns the data in the form of [time_soln1, timesoln2, ...] 
    specdata[species_name] = [spec_molefrac_soln1, spec_molefrac_soln2, ...].
    Time is returned in units of [seconds].
    """

    timeData = []
    distanceData = []
    specData = {spec: [] for spec in species}

    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            tokens = label.split('_')
            
            if tokens[0] == 'Time':
                tdata = np.array([float(r) for r in row[2:]], dtype=np.float64)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}.get(tunits, 1.0)
                timeData.append(tdata)
            
            elif tokens[0] == 'Distance':      
                ddata = np.array([float(r) for r in row[2:]], dtype=np.float64)
                dunits = row[1].strip()[1:-1].lower()
                ddata *= {'cm': 1.0, 'mm': 0.1, 'm': 100.}.get(dunits, 1.0)
                distanceData.append(ddata)

            elif label.startswith('Mole_fraction'):
                for spec in species:
                    if tokens[2] == spec:
                        specData[spec].append(np.array([float(r) for r in row[2:]], dtype=np.float64))

    if timeData:
        return timeData, specData
    elif distanceData:
        return distanceData, specData
    else:
        return


def getTotalMoles(ckcsvFile):
    """
    Return the time and total moles data from the given CKCSV file.  
    Returns the data in the form of [time_soln1, timesoln2, ...], [total_moles_soln1, total_moles_soln2, ...]
    Moles is returned in units of [moles].
    """
    
    tempData = []
    presData = []
    volData = []
    totalMolesData = []

    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            tokens = label.split('_')

            if tokens[0] == 'Temperature':
                Tdata = np.array([float(r) for r in row[2:]], dtype=np.float64)
                Tunits = row[1].strip()[1:-1].lower()
                try:
                    Tdata *= {'k': 1.0}[Tunits] 
                except KeyError:
                    print(f'Unknown units {Tunits} for temperature. Cannot obtain total moles.')
                tempData.append(Tdata)

            elif tokens[0] == 'Pressure': 
                Pdata = np.array([float(r) for r in row[2:]], dtype=np.float64)
                Punits = row[1].strip()[1:-1].lower()
                try:
                    Pdata *= {'bar': 0.1, 'mpa': 1.0, 'pa': 1e-6, 'atm': 0.101325}[Punits]
                except KeyError:
                    print(f'Unknown units {Punits} for pressure. Cannot obtain total moles.')
                presData.append(Pdata)
            
            elif tokens[0] == 'Volume':                
                Vdata = np.array([float(r) for r in row[2:]], dtype=np.float64)
                Vunits = row[1].strip()[1:-1].lower()
                try:
                    Vdata *= {'cm3': 1.0, 'm3': 1e6}[Vunits]
                except KeyError:
                    print(f'Unknown units {Vunits} for volume. Cannot obtain total moles.')
                volData.append(Vdata)

    R = 8.3145
    for i in range(len(tempData)):
        totalMolesData.append(presData[i] * volData[i] / (R * tempData[i]))

    return totalMolesData



################################################################################
def getIgnitionDelayOnset(csvfile,tol=1.0,secies =[],cond='onset'):
    """
    Return the ignition delay time from the given CKCSV file. This function
    uses (dP/dt)_max to find the ignition delay time. A ValueError is raised if
    this dP/dt value is below a certain threshold.
    
    Alternatively, provide a list of species that are maximized at tau.  
    If the list contains a single species label, e.g., ['OH'], the ignition delay will be given as the
    time at which the concentration of that species is maximized. If a list of species is given,
    then the multiplicative max of the concentrations of those species is returned. For instance, 
    if it is desired to calculate the ignition delay using max([CH][O]), then the species list 
    ['CH', 'O'] should be inputted.
    """
    
    tdata = None
    Pdata = None
    specdata = []

    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            if label.startswith('Time'):
                tdata = np.array([float(r) for r in row[2:]], dtype=np.float64)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}.get(tunits, 1.0)
            elif label.startswith('Pressure'):
                Pdata = np.array([float(r) for r in row[2:]], dtype=np.float64)
                Punits = row[1].strip()[1:-1].lower()
                Pdata *= {'dyne/cm2': 0.1, 'atm': 101325., 'Pa': 1.0, 'bar': 1e5, 'torr': 133.322368, 'mtorr': 0.133322368, 'psi': 6894.75729}.get(Punits, 1.0) * 1e-5
            for spec in species:
                if spec in label:
                    specdata.append(np.array([float(r) for r in row[2:]], dtype=np.float64))

    if tdata is None or Pdata is None:
        raise Exception('Unable to read time and/or pressure data from the given CKCSV file.')
    
    if species:
        if len(species) == 1:
            max_index = specdata[0].argmax()
            OHmetric = max(specdata[0]) / 2  # Reference value for specific mechanism
            mindata = abs(specdata[0][:max_index] - OHmetric)
            index = mindata.argmin()
            return tdata[index]
        else:
            multdata = np.ones(len(specdata[0]))
            ind1 = specdata[0].argmax()
            ind2 = specdata[1].argmax()
            for spec in specdata:  
                multdata *= spec
            print(multdata)
            print(tdata)
            print('time max for species 1')
            print(tdata[ind1])
            print('time max for species 2')
            print(tdata[ind2])
            index = multdata.argmax()
            print('time max for multiplied together')
            print(tdata[index])
            return tdata[index]

    else:
        #print(Pdata[1:],Pdata[:-1],tdata[1:],tdata[:-1])
        dPdt = (Pdata[1:] - Pdata[:-1]) / (tdata[1:] - tdata[:-1])
        dt = np.diff(tdata)
        dPdt = dPdt[np.isfinite(dPdt)]
        index = next(i for i, d in enumerate(dPdt) if d == max(dPdt))
        if dPdt[index] < tol:
            raise ValueError('No ignition occurred in the given simulation.')
       
        intercept = int(np.diff(dPdt).argmax())
        if cond == "dt-max":
            #print(tdata[index])
            return 0.5 * (tdata[index] + tdata[index + 1])
        else:
            return 0.5* (tdata[intercept-1] + tdata[intercept])

def getIgnitionDelay(ckcsvFile, tol=1.0, species=[],cond='onset'):
    """
    Return the ignition delay time from the given CKCSV file. This function
    uses (dP/dt)_max to find the ignition delay time. A ValueError is raised if
    this dP/dt value is below a certain threshold.
    
    Alternatively, provide a list of species that are maximized at tau.  
    If the list contains a single species label, e.g., ['OH'], the ignition delay will be given as the
    time at which the concentration of that species is maximized. If a list of species is given,
    then the multiplicative max of the concentrations of those species is returned. For instance, 
    if it is desired to calculate the ignition delay using max([CH][O]), then the species list 
    ['CH', 'O'] should be inputted.
    """
    
    tdata = None
    Pdata = None
    specdata = []

    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            if label.startswith('Time'):
                tdata = np.array([float(r) for r in row[2:]], dtype=np.float64)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}.get(tunits, 1.0)
            elif label.startswith('Pressure'):
                Pdata = np.array([float(r) for r in row[2:]], dtype=np.float64)
                Punits = row[1].strip()[1:-1].lower()
                Pdata *= {'dyne/cm2': 0.1, 'atm': 101325., 'Pa': 1.0, 'bar': 1e5, 'torr': 133.322368, 'mtorr': 0.133322368, 'psi': 6894.75729}.get(Punits, 1.0) * 1e-5
            for spec in species:
                if spec in label:
                    specdata.append(np.array([float(r) for r in row[2:]], dtype=np.float64))

    if tdata is None or Pdata is None:
        raise Exception('Unable to read time and/or pressure data from the given CKCSV file.')
    
    if species:
        if len(species) == 1:
            max_index = specdata[0].argmax()
            OHmetric = max(specdata[0]) / 2  # Reference value for specific mechanism
            mindata = abs(specdata[0][:max_index] - OHmetric)
            index = mindata.argmin()
            return tdata[index]
        else:
            multdata = np.ones(len(specdata[0]))
            ind1 = specdata[0].argmax()
            ind2 = specdata[1].argmax()
            for spec in specdata:  
                multdata *= spec
            print(multdata)
            print(tdata)
            print('time max for species 1')
            print(tdata[ind1])
            print('time max for species 2')
            print(tdata[ind2])
            index = multdata.argmax()
            print('time max for multiplied together')
            print(tdata[index])
            return tdata[index]

    else:
        #print(Pdata[1:],Pdata[:-1],tdata[1:],tdata[:-1])
        dPdt = (Pdata[1:] - Pdata[:-1]) / (tdata[1:] - tdata[:-1])
        dt = np.diff(tdata)
        dPdt = dPdt[np.isfinite(dPdt)]
        index = next(i for i, d in enumerate(dPdt) if d == max(dPdt))
        if dPdt[index] < tol:
            raise ValueError('No ignition occurred in the given simulation.')
       
        intercept = int(np.diff(dPdt).argmax())
        if cond == "dt-max":
            #print(tdata[index])
            return 0.5 * (tdata[index] + tdata[index + 1])
        else:
            return 0.5* (tdata[intercept-1] + tdata[intercept])

################################################################################

def runIgnitionReactionSensitivity(runChemkinJob, inputFile, dictionaryFile):
    """
    Supply a runChemkinJob python function which returns the ignition delay with a 
    Chemkin file input. This will run finite difference reaction sensitivities and 
    save them to a CSV file.
    """
    from rmgpy.chemkin import loadChemkinFile, saveChemkinFile
        
    speciesList, reactionList = loadChemkinFile(inputFile, dictionaryPath=dictionaryFile, readComments=False)
    num_reactions = len(reactionList)
    
    factor_high = 1.05
    factor_low = 1.0 / factor_high
    
    with open('ignition_rxn_sensitivity.csv', 'w', newline='') as csvfile:
        worksheet = csv.writer(csvfile)
        worksheet.writerow(['Index', 'Reaction', 'd[ln k]', 'tau_high', 'tau_low', 'd[ln tau]/d[ln k]'])
        
        logging.info('Running reaction sensitivity analysis using finite differences...')
        for index, reaction in enumerate(reactionList):
            rxn_index = index + 1
            rxn_string = reaction.toChemkin(kinetics=False)
            logging.info('At reaction {0} of {1}. {2}'.format(rxn_index, num_reactions, rxn_string))
            
            # Calculate tau_high
            reaction.kinetics.changeRate(factor_high)
            saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose=False)
            tau_high = runChemkinJob('chem_temp.inp')
            reaction.kinetics.changeRate(1.0 / factor_high)  # reset the kinetics
            
            # Calculate tau_low
            reaction.kinetics.changeRate(factor_low)
            saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose=False)
            tau_low = runChemkinJob('chem_temp.inp')
            reaction.kinetics.changeRate(1.0 / factor_low)   # reset the kinetics
            
            # Sensitivity calculation
            if tau_high != 0 and tau_low != 0:
                sens = np.log(tau_high / tau_low) / np.log(factor_high / factor_low)
            else:
                sens = 0
            
            worksheet.writerow([rxn_index, rxn_string, np.log(factor_high / factor_low), tau_high, tau_low, sens])

        
################################################################################

def runIgnitionThermoSensitivity(runChemkinJob, inputFile, dictionaryFile):
    """
    Supply a runChemkinJob python function which returns the ignition delay with a 
    Chemkin file input. This will run finite difference sensitivities to enthalpies and 
    save them to a CSV file.
    """
    
    from rmgpy.chemkin import loadChemkinFile, saveChemkinFile, getSpeciesIdentifier
    from rmgpy.quantity import Quantity
    
    speciesList, reactionList = loadChemkinFile(inputFile, dictionaryPath=dictionaryFile, readComments=False)
    num_species = len(speciesList)
    
    deltaH = Quantity(0.5, 'kcal/mol').value_si
    
    with open('ignition_thermo_sensitivity.csv', 'w', newline='') as csvfile:
        worksheet = csv.writer(csvfile)
        worksheet.writerow(['Species', 'd[del H] (kcal/mol)', 'tau_high', 'tau_low', 'd[ln tau]/d[del H]'])
        
        logging.info('Running thermo sensitivity analysis using finite differences...')
        for index, species in enumerate(speciesList):
            species_index = index + 1
            species_string = getSpeciesIdentifier(species)
            logging.info(f'At species {species_index} of {num_species}. {species_string}')
            
            # Increase enthalpy
            species.thermo.changeBaseEnthalpy(deltaH)
            saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose=False)
            tau_high = runChemkinJob('chem_temp.inp')
            species.thermo.changeBaseEnthalpy(-deltaH)  # reset the thermo
            
            # Decrease enthalpy
            species.thermo.changeBaseEnthalpy(-deltaH)
            saveChemkinFile('chem_temp.inp', speciesList, reactionList, verbose=False)
            tau_low = runChemkinJob('chem_temp.inp')
            species.thermo.changeBaseEnthalpy(deltaH)   # reset the thermo
            
            # Sensitivity calculation
            if tau_high != 0 and tau_low != 0:
                sens = np.log(tau_high / tau_low) / (2 * deltaH)
            else:
                sens = 0
                
            worksheet.writerow([species_string, '1', tau_high, tau_low, sens])

def getIgnitionDelayOH(ckcsvFile, tol=1.0):
    """
    Return the ignition delay time from the given CKCSV file. This function
    uses (OH)_max to find the ignition delay time. 
    """
    
    tdata = None
    OHdata = None
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            if label.startswith('Time'):
                tdata = np.array([float(r) for r in row[2:]], dtype=float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
            elif label.startswith('Mole_fraction_OH'):
                OHdata = np.array([float(r) for r in row[2:]], dtype=float)
    
    if tdata is None or OHdata is None:
        raise Exception('Unable to read time and/or OH data from the given CKCSV file.')
    
    OHdata = OHdata[np.isfinite(OHdata)] 
    index = OHdata.argmax()

    if index == len(OHdata) - 1:
        raise ValueError('No ignition occurred in the given simulation.')

    return tdata[index]


################################################################################



def getPeakOQOOHTime(ckcsvFile, spc, tol=1.0):
    """
    Return the ignition delay time from the given CKCSV file. This function
    uses (OQOOH)_max to find the ignition delay time based on the peak mole fraction 
    of the specified species (OQOOH).
    """
    
    import numpy as np
    
    tdata = None
    OQOOHdata = None
    spc_label = 'Mole_fraction_' + spc
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            if label.startswith('Time'):
                tdata = np.array([float(r) for r in row[2:]], dtype=float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
            elif label.startswith(spc_label):
                OQOOHdata = np.array([float(r) for r in row[2:]], dtype=float)
    
    if tdata is None or OQOOHdata is None:
        raise Exception('Unable to read time and/or OQOOH data from the given CKCSV file.')
    
    OQOOHdata = OQOOHdata[np.isfinite(OQOOHdata)]
    index = OQOOHdata.argmax()

    if index == len(OQOOHdata) - 1:
        raise ValueError('No OQOOH peak found in the given simulation.')

    return tdata[index]


################################################################################

def getIgnitionDelayStage1(ckcsvFile, stepsize=1500, tol=1.0):
    """
    Return the first-stage ignition delay time from the given CKCSV file.
    This function uses (d²T/dt²) inflection to find the first stage ignition delay time.
    A ValueError is raised if the dT/dt inflection does not meet the threshold.
    """
    import numpy as np
    from scipy import interpolate

    tdata = None
    Tdata = None
    
    with open(ckcsvFile, 'r') as stream:
        reader = csv.reader(stream)
        for row in reader:
            label = row[0].strip()
            if label.startswith('Time'):
                tdata = np.array([float(r) for r in row[2:]], dtype=float)
                tunits = row[1].strip()[1:-1].lower()
                tdata *= {'sec': 1.0, 'min': 60., 'hr': 3600., 'msec': 1e-3, 'microsec': 1e-6}[tunits]
            elif label.startswith('Temperature'):
                Tdata = np.array([float(r) for r in row[2:]], dtype=float)
    
    if tdata is None or Tdata is None:
        raise Exception('Unable to read time and/or temperature data from the given CKCSV file.')
    
    # Interpolating the temperature data over a defined range
    temperature_interp = interpolate.interp1d(tdata, Tdata, kind='cubic')
    xs = np.linspace(0.0, tdata[-1], stepsize)
    ys = temperature_interp(xs)
    
    # Calculating first and second derivatives with respect to time
    dT_dt = np.diff(ys) / np.diff(xs)
    xs1 = xs[1:]
    dT2_dt2 = np.diff(dT_dt) / np.diff(xs1)
    xs2 = xs1[1:]
    
    # Finding the first significant inflection where d²T/dt² < -tol
    try:
        time_index_Tinflection = next(i for i, val in enumerate(dT2_dt2) if val < -tol)
    except StopIteration:
        raise ValueError('No T kink found in the data.')
    
    Time_inflection = xs2[time_index_Tinflection]

    return Time_inflection

################################################################################
