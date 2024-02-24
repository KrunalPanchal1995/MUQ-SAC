# MUQ-SAC
The suit of code is part of project undertaken for kinetic mechanism optimization framework development. The code includes the temperature-dependent Modified Uncertainty Quantification (MUQ) method and Sampling of Arrhenius Curves (SAC) technique.

Step 1

 - Converting mechanism files into yaml: (write code in terminal)
     
     ck2yaml --input=chem.inp --thermo=therm.dat --transport=tran.dat --permissive
     (eg:ck2yaml --input=MB_MB2D.inp --thermo=MB+MB2D_therm.dat --transport=Mb+MB2d_trans.dat --permissive)
     
Step 2

 - Generating xml file for reactions considered in sensitivity analysis
 
Step 3
 
 - Generating target input file

