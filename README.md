# MUQ-SAC: Methof of Uncertainty Quantification and Sampling of ARRHENIUS CURVES
The suit of code is part of project undertaken for kinetic mechanism optimization framework development. The code includes the temperature-dependent Modified Uncertainty Quantification (MUQ) method and Sampling of Arrhenius Curves (SAC) technique.

Step 1

 - Converting mechanism files into yaml: (write code in terminal)
     
     ck2yaml --input=chem.inp --thermo=therm.dat --transport=tran.dat --permissive
     (eg:ck2yaml --input=MB_MB2D.inp --thermo=MB+MB2D_therm.dat --transport=Mb+MB2d_trans.dat --permissive)
     
Step 2

 - Generating xml file for reactions considered in sensitivity analysis
 
Step 3
 
 - Generating target input file [target file, Addendum file]
 - Collect the mechanism file [mech file, thermo file, transport file]

#For Optimization
Step 4
- run in terminal
	python3.9 /home/krithika/Desktop/MUQ-SAC/MUQ-SAC/run_script.py [target_file]



#For running simple simulations
Step 4
- run in terminal
	python3.9 /home/krithika/Desktop/MUQ-SAC/MUQ-SAC/run_nominal_sim.py [target_file]


#For running sensitivity analysis
Step 4
- run in terminal
	python3.9 /home/krithika/Desktop/MUQ-SAC/SENSITIVITY_ANALYSIS/sens.py [target_file]
 

