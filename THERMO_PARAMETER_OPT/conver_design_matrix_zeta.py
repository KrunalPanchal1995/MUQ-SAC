
import numpy as np
import pickle
import os
import sys
import pandas as pd
import numpy as np

# Uncertainty data (nominal coefficients, covariance matrix, and zeta_max)
# The covariance matrices are stored as numpy arrays for matrix operations.

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: AR - Low Temperature Range
#--------------------------------------------------------------------------------

AR_Low_nominal_param = [
     2.500000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00    ]

AR_Low_cov_matrix = np.array([
    [ 2.139158220517e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-4.022124051830e-05,  2.845717273878e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-8.890603798845e-08, -9.974455166573e-08, -6.903561866731e-09,  0.000000000000e+00,  0.000000000000e+00],
    [ 1.192900331117e-12,  2.636679886098e-12,  5.883533113896e-12,  8.341224950427e-12,  0.000000000000e+00],
    [ 1.371794841382e-15,  1.512386609301e-15,  1.509735912439e-15,  1.621523761220e-15,  1.605197991945e-15],
])

AR_Low_max_zeta = [
     2.231625019310e+00,  2.120291290484e+00,  4.225985980384e-04,  6.347857042679e-02, 
     1.729694231522e-02
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: AR - High Temperature Range
#--------------------------------------------------------------------------------

AR_High_nominal_param = [
     2.500000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,    0.000000000000e+00    ]

AR_High_cov_matrix = np.array([
    [-2.481478181167e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 8.990734069841e-05,  9.964147556348e-05,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-4.801392766627e-08, -3.191490927450e-08, -7.899689954706e-09,  0.000000000000e+00,  0.000000000000e+00],
    [ 9.210544867977e-12, -3.255074569015e-13,  4.536715791933e-12,  5.617678736556e-12,  0.000000000000e+00],
    [-3.738698626471e-16,  3.234911590340e-16, -5.226571390029e-16, -9.692338051238e-16,  1.308651009837e-16],
])

AR_High_max_zeta = [
    -3.010865949532e+00,  5.451194286620e-01, -1.202506357046e-01, -4.005970534850e-01, 
     6.182830375093e+00
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: H2 - Low Temperature Range
#--------------------------------------------------------------------------------

H2_Low_nominal_param = [
     2.344331120000e+00,  7.980520750000e-03, -1.947815100000e-05,  2.015720940000e-08, 
    -7.376117610000e-12    ]

H2_Low_cov_matrix = np.array([
    [-3.125890217219e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 1.679651679878e-04,  6.134776628442e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 1.724788777132e-07, -3.102364670871e-07,  3.573419545815e-08,  0.000000000000e+00,  0.000000000000e+00],
    [ 1.588281597909e-11, -6.971389790468e-12,  7.737441419145e-12,  8.269571054768e-12,  0.000000000000e+00],
    [ 1.921043894675e-15,  1.184610087582e-15,  1.593786312437e-15,  1.615908276981e-15,  1.603590521149e-15],
])

H2_Low_max_zeta = [
    -2.708220785513e+00,  9.757493197554e-01,  1.542217230544e+01,  4.846715265558e+00, 
     9.663908359478e-01
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: H2 - High Temperature Range
#--------------------------------------------------------------------------------

H2_High_nominal_param = [
     2.932865750000e+00,  8.266080260000e-04, -1.464023640000e-07,  1.541004140000e-11, 
    -6.888048000000e-16    ]

H2_High_cov_matrix = np.array([
    [-3.862200372162e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 3.388671696972e-04,  5.702836849286e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 5.210755824255e-08, -3.143475427760e-07, -1.034137280709e-08,  0.000000000000e+00,  0.000000000000e+00],
    [-3.762460590986e-11,  4.691961735121e-11,  1.786020496248e-12,  3.397538434832e-12,  0.000000000000e+00],
    [ 3.473151397913e-15, -1.948190044860e-15,  1.211349901507e-16, -4.047611557763e-16,  9.710396212918e-18],
])

H2_High_max_zeta = [
    -2.137859243547e-01,  2.437493134269e+00, -9.760134974452e+00,  1.447317148123e+01, 
     2.000130998144e-02
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: O2 - Low Temperature Range
#--------------------------------------------------------------------------------

O2_Low_nominal_param = [
     3.782456360000e+00, -2.996734160000e-03,  9.847302010000e-06, -9.681295090000e-09, 
     3.243728370000e-12
]

O2_Low_cov_matrix = np.array([
    [-3.115697959512e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 1.872271483717e-04,  5.850900335417e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 3.824637440004e-08, -2.382785270033e-07,  1.594792511335e-08,  0.000000000000e+00,  0.000000000000e+00],
    [ 7.708851030716e-12, -8.784905791949e-12,  6.668830296348e-12,  7.602691349972e-12,  0.000000000000e+00],
    [ 1.576371989804e-15,  8.225356039927e-16,  1.539427910283e-15,  1.583481665389e-15,  1.597181130591e-15],
])

O2_Low_max_zeta = [
    -2.377364852194e+00,  1.678652969218e+00,  1.042289509906e+01,  5.245236059017e+00, 
     1.097700520128e+00
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: O2 - High Temperature Range
#--------------------------------------------------------------------------------

O2_High_nominal_param = [
     3.660960650000e+00,  6.563658110000e-04, -1.411496270000e-07,  2.057979350000e-11, 
    -1.299134360000e-15
]

O2_High_cov_matrix = np.array([
    [ 3.059487939289e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 5.086837416012e-05, -5.258999012284e-05,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-1.096368684554e-08, -8.024944978513e-09, -1.588272054630e-08,  0.000000000000e+00,  0.000000000000e+00],
    [-2.127789730745e-12,  2.593830047143e-13,  1.779992056277e-12, -9.458756720668e-13,  0.000000000000e+00],
    [ 3.474159126016e-16,  1.997668119585e-16,  4.729925032192e-17,  2.190303858533e-16,  1.930907677500e-16],
])

O2_High_max_zeta = [
     2.781216107142e+00, -1.107994549200e+00, -3.620105137844e-01,  6.443208616025e-01, 
     2.040212606654e+00
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: H2O - Low Temperature Range
#--------------------------------------------------------------------------------

H2O_Low_nominal_param = [
     4.198635200000e+00, -2.036401700000e-03,  6.520341600000e-06, -5.487926900000e-09, 
     1.771968000000e-12
]

H2O_Low_cov_matrix = np.array([
    [-3.561180303794e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 1.757248066403e-04,  5.781386946531e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 4.901067268364e-08, -1.804342052877e-07,  2.375049249320e-08,  0.000000000000e+00,  0.000000000000e+00],
    [ 8.556165219868e-12, -5.617940534252e-12,  7.026395499361e-12,  7.731576377324e-12,  0.000000000000e+00],
    [ 1.622108747670e-15,  9.595512890601e-16,  1.554193039153e-15,  1.588580343724e-15,  1.598034305590e-15],
])

H2O_Low_max_zeta = [
    -2.561124310072e+00,  1.345034684374e+00,  1.118987977731e+01,  4.643104687789e+00, 
     9.810365126131e-01
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: H2O - High Temperature Range
#--------------------------------------------------------------------------------

H2O_High_nominal_param = [
     2.677038900000e+00,  2.973181600000e-03, -7.737688900000e-07,  9.443351400000e-11, 
    -4.268999100000e-15
]

H2O_High_cov_matrix = np.array([
    [-3.461779300156e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 6.688767307947e-05, -3.396528612694e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-4.130439383135e-08,  7.475386800694e-08,  1.765297819484e-08,  0.000000000000e+00,  0.000000000000e+00],
    [-5.082363648126e-12,  7.759550535635e-12,  2.631408259955e-12,  7.057632550209e-12,  0.000000000000e+00],
    [ 2.122194922500e-15, -1.202700308904e-15, -8.546178799364e-16, -1.123516269122e-15,  1.655697365662e-16],
])

H2O_High_max_zeta = [
    -2.402152766068e+00, -1.834343987646e+00,  4.351584696016e-01, -1.992736204390e+00, 
     1.895569387143e+01
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: C - Low Temperature Range
#--------------------------------------------------------------------------------

C_Low_nominal_param = [
     2.554239500000e+00, -3.215377200000e-04,  7.337922300000e-07, -7.322348700000e-10, 
     2.665214400000e-13
]

C_Low_cov_matrix = np.array([
    [ 2.152248591931e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-4.382695425675e-05,  2.841320258874e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-8.546660796831e-08, -9.972817885862e-08, -6.788268130096e-09,  0.000000000000e+00,  0.000000000000e+00],
    [ 1.338724723191e-12,  2.484853593735e-12,  5.887071715625e-12,  8.319308601288e-12,  0.000000000000e+00],
    [ 1.374644577310e-15,  1.498882427205e-15,  1.509835297822e-15,  1.620505568774e-15,  1.604996861748e-15],
])

C_Low_max_zeta = [
     2.238291529958e+00,  2.111816022313e+00,  1.268940160609e-03,  6.324411823858e-02, 
     1.728021045140e-02
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: C - High Temperature Range
#--------------------------------------------------------------------------------

C_High_nominal_param = [
     2.605583000000e+00, -1.959343400000e-04,  1.067372200000e-07, -1.642394000000e-11, 
     8.187058000000e-16
]

C_High_cov_matrix = np.array([
    [-2.106303286833e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 2.544452795121e-06, -3.157378358891e-05,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 4.187553732769e-09,  5.488572882083e-09, -9.998233996811e-09,  0.000000000000e+00,  0.000000000000e+00],
    [-3.106138493547e-12,  2.740542713679e-13,  5.326417424696e-12,  1.051220671900e-12,  0.000000000000e+00],
    [ 4.663824674306e-16, -5.735239243455e-17, -5.832345198145e-16, -1.786793000884e-16,  1.834208543534e-16],
])

C_High_max_zeta = [
    -3.039615149852e+00,  3.435756785966e-02, -5.639802234294e-01, -6.651531124854e-01, 
     1.298964578525e+00
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: CO - Low Temperature Range
#--------------------------------------------------------------------------------

CO_Low_nominal_param = [
     3.579533500000e+00, -6.103536900000e-04,  1.016814300000e-06,  9.070058600000e-10, 
    -9.044244900000e-13
]

CO_Low_cov_matrix = np.array([
    [-3.357020915409e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 2.754567892646e-04,  5.877394209526e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 3.181297971732e-08, -2.489473424743e-07,  1.495013101790e-08,  0.000000000000e+00,  0.000000000000e+00],
    [ 7.087308858117e-12, -9.054622204699e-12,  6.633592830719e-12,  7.593357083680e-12,  0.000000000000e+00],
    [ 1.542005925677e-15,  8.226603993461e-16,  1.538561219648e-15,  1.583267230391e-15,  1.597162267436e-15],
])

CO_Low_max_zeta = [
    -2.376595977369e+00,  1.620080984838e+00,  1.507029796626e+01,  7.857144593784e+00, 
     1.643663315529e+00
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: CO - High Temperature Range
#--------------------------------------------------------------------------------

CO_High_nominal_param = [
     3.048485900000e+00,  1.351728100000e-03, -4.857940500000e-07,  7.885364400000e-11, 
    -4.698074600000e-15
]

CO_High_cov_matrix = np.array([
    [-4.315167023834e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 3.575192106823e-04, -4.908542411782e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-8.183069983357e-08,  2.484507268064e-07,  7.511632505507e-08,  0.000000000000e+00,  0.000000000000e+00],
    [-2.938805451686e-12, -4.886188904716e-11, -1.796459018362e-11,  3.746520447765e-12,  0.000000000000e+00],
    [ 1.264966451441e-15,  3.831500090996e-15,  8.822515045011e-16, -4.315035129214e-16, -3.400833660292e-17],
])

CO_High_max_zeta = [
    -1.523142570646e-01, -3.303107279292e+00, -5.987140672256e-01,  6.470648245361e+00, 
    -1.305385169942e+00
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: OH* - Low Temperature Range
#--------------------------------------------------------------------------------
# Note: Python variable name uses 'OHstar' for 'OH*'

OHstar_Low_nominal_param = [
     3.460844280000e+00,  5.018721720000e-04, -2.002544740000e-06,  3.189019840000e-09, 
    -1.354518380000e-12
]

OHstar_Low_cov_matrix = np.array([
    [-3.251477935931e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 2.108830054753e-04,  5.846724448517e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 1.189405876295e-07, -2.663519705465e-07,  3.388563977785e-08,  0.000000000000e+00,  0.000000000000e+00],
    [ 1.248248563123e-11, -6.150791169772e-12,  7.673578145506e-12,  8.096437057572e-12,  0.000000000000e+00],
    [ 1.770378553226e-15,  1.143558786908e-15,  1.591263863060e-15,  1.607710406496e-15,  1.601957281939e-15],
])

OHstar_Low_max_zeta = [
    -2.628488245231e+00, 1.103131375376e+00, 1.417044402381e+01, 4.574125732498e+00, 9.258389746243e-01
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: OH* - High Temperature Range
# (Note: Python variable name uses 'OHstar' for 'OH*')
#--------------------------------------------------------------------------------

OHstar_High_nominal_param = [
     2.755829200000e+00,  1.398487560000e-03, -4.194284930000e-07,  6.334532820000e-11, 
    -3.560422180000e-15
]

OHstar_High_cov_matrix = np.array([
    [ 2.555623434587e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 5.051012144146e-05,  1.108231139799e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 7.129850261146e-09, -4.710163492070e-08, -4.432269775558e-08,  0.000000000000e+00,  0.000000000000e+00],
    [-4.140559184844e-12, -4.071812494329e-12,  1.407136982620e-11,  3.703597713101e-12,  0.000000000000e+00],
    [-1.137142096217e-16,  1.220777718175e-15, -1.038707448883e-15, -6.088115146579e-16,  2.143654475433e-16],
])

OHstar_High_max_zeta = [
     3.157864042835e+00, -8.722102870577e-02,  3.973498124700e-01, -9.513151859393e-01, 
     7.316740034899e+00
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: O - Low Temperature Range
#--------------------------------------------------------------------------------

O_Low_nominal_param = [
     3.168267100000e+00, -3.279318840000e-03,  6.643063960000e-06, -6.128066240000e-09, 
     2.112659710000e-12
]

O_Low_cov_matrix = np.array([
    [ 2.425778881168e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-1.161414382712e-04,  2.817612149263e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-2.284627608996e-08, -1.033306821696e-07,  6.523364802696e-09,  0.000000000000e+00,  0.000000000000e+00],
    [ 3.795850332885e-12, -1.487441167948e-12,  6.237926760162e-12,  7.445587331996e-12,  0.000000000000e+00],
    [ 1.387501699950e-15,  1.126626516400e-15,  1.518983212828e-15,  1.576441934018e-15,  1.595887343209e-15],
])

O_Low_max_zeta = [
     2.362987675108e+00,  1.947079219119e+00,  1.007871546898e-01,  6.181522956512e-02, 
     1.819281603588e-02
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: O - High Temperature Range
#--------------------------------------------------------------------------------

O_High_nominal_param = [
     2.543636970000e+00, -2.731624860000e-05, -4.190295200000e-09,  4.954818450000e-12, 
    -4.795536940000e-16
]

O_High_cov_matrix = np.array([
    [ 2.624326349901e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-1.707995674924e-04,  3.665165436272e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-5.138994933896e-08, -2.007263600402e-07,  2.855433236066e-08,  0.000000000000e+00,  0.000000000000e+00],
    [ 2.889266734101e-11,  3.398828857804e-11, -8.154631340832e-12,  6.904203456003e-12,  0.000000000000e+00],
    [-2.885355962038e-15, -1.977845936370e-15,  5.034068792405e-16, -1.202634723274e-15, -4.067446901861e-17],
])

O_High_max_zeta = [
    -1.223814582581e+00,  3.893597258075e+00, -8.271892059807e+00,  9.110086982129e+00, 
    -1.000133801928e+01
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: HCO - Low Temperature Range
#--------------------------------------------------------------------------------

HCO_Low_nominal_param = [
     4.237546100000e+00, -3.320752570000e-03,  1.400302640000e-05, -1.342399950000e-08, 
     4.374162080000e-12
]

HCO_Low_cov_matrix = np.array([
    [-2.697559915232e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],                
    [-2.605694645926e-04,  9.903892638336e-05,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 4.759252532127e-08,  2.322184726505e-08,  3.285990322897e-08,  0.000000000000e+00,  0.000000000000e+00],
    [ 9.548344337882e-12,  6.886863132766e-12,  7.629623815474e-12,  7.916770773049e-12,  0.000000000000e+00],
    [ 1.703572116513e-15,  1.542085299300e-15,  1.582782596920e-15,  1.596350028709e-15,  1.599360902737e-15],
])

HCO_Low_max_zeta = [
    -2.943314107684e+00,  5.858933757016e-01,  1.503765747609e-01,  3.567368263533e-02, 
     1.368372424301e-02
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: HCO - High Temperature Range
#--------------------------------------------------------------------------------

HCO_High_nominal_param = [
     3.920015420000e+00,  2.522793240000e-03, -6.710041640000e-07,  1.056159480000e-10, 
    -7.437982610000e-15
]

HCO_High_cov_matrix = np.array([
    [-3.453881188140e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-1.648683021937e-04,  1.205593280795e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 3.486312749030e-08,  5.435244244457e-09, -5.668188168546e-09,  0.000000000000e+00,  0.000000000000e+00],
    [ 3.597226961796e-12,  1.922703761432e-12, -5.984628323456e-13,  6.574492378742e-12,  0.000000000000e+00],
    [-6.666677171186e-16, -6.942617101103e-16,  3.167200055184e-16, -9.397966269010e-16,  3.826735497302e-16],
])

HCO_High_max_zeta = [
    -3.029838622596e+00,  1.571747774036e-01, -1.580207069411e-01,  3.357432106778e+00, 
     1.677432312071e+00
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: H - Low Temperature Range
# (Note: This data block is identical to the AR Low-T data provided previously.)
#--------------------------------------------------------------------------------

H_Low_nominal_param = [
     2.500000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00, 
     0.000000000000e+00
]

H_Low_cov_matrix = np.array([
    [ 2.139158220517e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-4.022124051830e-05,  2.845717273878e-04,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-8.890603798845e-08, -9.974455166573e-08, -6.903561866731e-09,  0.000000000000e+00,  0.000000000000e+00],
    [ 1.192900331117e-12,  2.636679886098e-12,  5.883533113896e-12,  8.341224950427e-12,  0.000000000000e+00],
    [ 1.371794841382e-15,  1.512386609301e-15,  1.509735912439e-15,  1.621523761220e-15,  1.605197991945e-15],
])

H_Low_max_zeta = [
     2.231625019310e+00,  2.120291290484e+00,  4.225985980384e-04,  6.347857042679e-02, 
     1.729694231522e-02
]

#--------------------------------------------------------------------------------
# UNCERTAINTY DATA FOR SPECIES: H - High Temperature Range
# (Note: This data block is identical to the AR High-T data provided previously.)
#--------------------------------------------------------------------------------

H_High_nominal_param = [
     2.500000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00, 
     0.000000000000e+00
]

H_High_cov_matrix = np.array([
    [-2.481478181167e-01,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [ 8.990734069841e-05,  9.964147556348e-05,  0.000000000000e+00,  0.000000000000e+00,  0.000000000000e+00],
    [-4.801392766627e-08, -3.191490927450e-08, -7.899689954706e-09,  0.000000000000e+00,  0.000000000000e+00],
    [ 9.210544867977e-12, -3.255074569015e-13,  4.536715791933e-12,  5.617678736556e-12,  0.000000000000e+00],
    [-3.738698626471e-16,  3.234911590340e-16, -5.226571390029e-16, -9.692338051238e-16,  1.308651009837e-16],
])

H_High_max_zeta = [
    -3.010865949532e+00,  5.451194286620e-01, -1.202506357046e-01, -4.005970534850e-01, 
     6.182830375093e+00
]


import numpy as np
from typing import Union, List, Tuple

import numpy as np

import numpy as np
from typing import Union, List, Tuple

# NOTE: This function assumes that the 40 variables (e.g., AR_Low_nominal_param  , AR_Low_cov_matrix, etc.) 
# have already been defined globally as NumPy arrays (5-element vector and 5x5 matrix, respectively)
# in the script where this function is executed.

def process_dm_zeta(input_vector):
    """
    Takes a 100-element perturbation vector, validates its length, slices it 
    into 20 blocks, and calculates the 20 modified parameter vectors.
    """
    
    # Ensure input is a numpy array (or can be converted) for robust slicing
    dm_zeta_vector = np.array(input_vector).flatten()
    
    # --- Step 1: Validation and Slicing (As previously defined) ---
    
    if len(dm_zeta_vector) != 100:
        raise AssertionError(f"Input vector must have exactly 100 elements, but found {len(dm_zeta_vector)}.")
    
    # Slicing and Assignment
    
    # Species 1: AR (Indices 0 to 9)
    AR_Low_dm_zeta = dm_zeta_vector[0:5]
    AR_High_dm_zeta = dm_zeta_vector[5:10]

    # Species 2: H2 (Indices 10 to 19)
    H2_Low_dm_zeta = dm_zeta_vector[10:15]
    H2_High_dm_zeta = dm_zeta_vector[15:20]

    # Species 3: O (Indices 20 to 29)
    O_Low_dm_zeta = dm_zeta_vector[20:25]
    O_High_dm_zeta = dm_zeta_vector[25:30]

    # Species 4: O2 (Indices 30 to 39)
    O2_Low_dm_zeta = dm_zeta_vector[30:35]
    O2_High_dm_zeta = dm_zeta_vector[35:40]

    # Species 5: H2O (Indices 40 to 49)
    H2O_Low_dm_zeta = dm_zeta_vector[40:45]
    H2O_High_dm_zeta = dm_zeta_vector[45:50]

    # Species 6: CO (Indices 50 to 59)
    CO_Low_dm_zeta = dm_zeta_vector[50:55]
    CO_High_dm_zeta = dm_zeta_vector[55:60]

    # Species 7: C (Indices 60 to 69)
    C_Low_dm_zeta = dm_zeta_vector[60:65]
    C_High_dm_zeta = dm_zeta_vector[65:70]

    # Species 8: HCO (Indices 70 to 79)
    HCO_Low_dm_zeta = dm_zeta_vector[70:75]
    HCO_High_dm_zeta = dm_zeta_vector[75:80]

    # Species 9: OH* (Indices 80 to 89) (Using OHstar)
    OHstar_Low_dm_zeta = dm_zeta_vector[80:85]
    OHstar_High_dm_zeta = dm_zeta_vector[85:90]

    # Species 10: H (Indices 90 to 99)
    H_Low_dm_zeta = dm_zeta_vector[90:95]
    H_High_dm_zeta = dm_zeta_vector[95:100]
    
    # --- Step 2: Calculate Modified Parameters Explicitly ---
    '''
    # Formula: modified_X = X_p0 + X_Cov * X_dm_zeta (where * is matrix dot product)
    
    # Species 1: AR
    modified_AR_Low_dm_zeta = AR_Low_nominal_param   + np.dot(AR_Low_cov_matrix, AR_Low_dm_zeta)
    modified_AR_High_dm_zeta = AR_High_nominal_param + np.dot(AR_High_cov_matrix, AR_High_dm_zeta)

    # Species 2: H2
    modified_H2_Low_dm_zeta = H2_Low_nominal_param   + np.dot(H2_Low_cov_matrix, H2_Low_dm_zeta)
    modified_H2_High_dm_zeta = H2_High_nominal_param + np.dot(H2_High_cov_matrix, H2_High_dm_zeta)

    # Species 3: O
    modified_O_Low_dm_zeta = O_Low_nominal_param   + np.dot(O_Low_cov_matrix, O_Low_dm_zeta)
    modified_O_High_dm_zeta = O_High_nominal_param + np.dot(O_High_cov_matrix, O_High_dm_zeta)

    # Species 4: O2
    modified_O2_Low_dm_zeta = O2_Low_nominal_param   + np.dot(O2_Low_cov_matrix, O2_Low_dm_zeta)
    modified_O2_High_dm_zeta = O2_High_nominal_param + np.dot(O2_High_cov_matrix, O2_High_dm_zeta)

    # Species 5: H2O
    modified_H2O_Low_dm_zeta = H2O_Low_nominal_param   + np.dot(H2O_Low_cov_matrix, H2O_Low_dm_zeta)
    modified_H2O_High_dm_zeta = H2O_High_nominal_param + np.dot(H2O_High_cov_matrix, H2O_High_dm_zeta)

    # Species 6: CO
    modified_CO_Low_dm_zeta = CO_Low_nominal_param   + np.dot(CO_Low_cov_matrix, CO_Low_dm_zeta)
    modified_CO_High_dm_zeta = CO_High_nominal_param + np.dot(CO_High_cov_matrix, CO_High_dm_zeta)

    # Species 7: C
    modified_C_Low_dm_zeta = C_Low_nominal_param   + np.dot(C_Low_cov_matrix, C_Low_dm_zeta)
    modified_C_High_dm_zeta = C_High_nominal_param + np.dot(C_High_cov_matrix, C_High_dm_zeta)

    # Species 8: HCO
    modified_HCO_Low_dm_zeta = HCO_Low_nominal_param   + np.dot(HCO_Low_cov_matrix, HCO_Low_dm_zeta)
    modified_HCO_High_dm_zeta = HCO_High_nominal_param + np.dot(HCO_High_cov_matrix, HCO_High_dm_zeta)

    # Species 9: OH*
    modified_OHstar_Low_dm_zeta = OHstar_Low_nominal_param   + np.dot(OHstar_Low_cov_matrix, OHstar_Low_dm_zeta)
    modified_OHstar_High_dm_zeta = OHstar_High_nominal_param + np.dot(OHstar_High_cov_matrix, OHstar_High_dm_zeta)

    # Species 10: H
    modified_H_Low_dm_zeta = H_Low_nominal_param   + np.dot(H_Low_cov_matrix, H_Low_dm_zeta)
    modified_H_High_dm_zeta = H_High_nominal_param + np.dot(H_High_cov_matrix, H_High_dm_zeta)
    
    
    # --- Step 3: Concatenate the 20 modified variables into a single 100-element vector ---
    modified_dm_zeta_vector = np.concatenate([
        modified_AR_Low_dm_zeta, modified_AR_High_dm_zeta,
        modified_H2_Low_dm_zeta, modified_H2_High_dm_zeta,
        modified_O_Low_dm_zeta, modified_O_High_dm_zeta,
        modified_O2_Low_dm_zeta, modified_O2_High_dm_zeta,
        modified_H2O_Low_dm_zeta, modified_H2O_High_dm_zeta,
        modified_CO_Low_dm_zeta, modified_CO_High_dm_zeta,
        modified_C_Low_dm_zeta, modified_C_High_dm_zeta,
        modified_HCO_Low_dm_zeta, modified_HCO_High_dm_zeta,
        modified_OHstar_Low_dm_zeta, modified_OHstar_High_dm_zeta,
        modified_H_Low_dm_zeta, modified_H_High_dm_zeta
    ])
    
    # Return the final 100-entry modified vector
    return modified_dm_zeta_vector    '''
    
    # CONVERT x_i ( ROW OF DESIGN MATRIX) to z_i )
    theta_low = [ 300, 500, 700, 800, 1000]
    theta_high= [1000, 1500, 2000, 3000, 5000]
    
    
    # Species 1: AR
    # NOTE: The AR_Low/High_max_zeta argument is added to match the 5-argument function signature.
    modified_AR_Low_dm_zeta = calculate_z_vector(theta_low, AR_Low_nominal_param, AR_Low_cov_matrix, AR_Low_max_zeta, AR_Low_dm_zeta)
    modified_AR_High_dm_zeta = calculate_z_vector(theta_high, AR_High_nominal_param, AR_High_cov_matrix, AR_High_max_zeta, AR_High_dm_zeta)

    # Species 2: H2
    modified_H2_Low_dm_zeta = calculate_z_vector(theta_low, H2_Low_nominal_param, H2_Low_cov_matrix, H2_Low_max_zeta, H2_Low_dm_zeta)
    modified_H2_High_dm_zeta = calculate_z_vector(theta_high, H2_High_nominal_param, H2_High_cov_matrix, H2_High_max_zeta, H2_High_dm_zeta)

    # Species 3: O (Now using the inverse formula)
    modified_O_Low_dm_zeta = calculate_z_vector(theta_low, O_Low_nominal_param, O_Low_cov_matrix, O_Low_max_zeta, O_Low_dm_zeta)
    modified_O_High_dm_zeta = calculate_z_vector(theta_high, O_High_nominal_param, O_High_cov_matrix, O_High_max_zeta, O_High_dm_zeta)

    # Species 4: O2
    modified_O2_Low_dm_zeta = calculate_z_vector(theta_low, O2_Low_nominal_param, O2_Low_cov_matrix, O2_Low_max_zeta, O2_Low_dm_zeta)
    modified_O2_High_dm_zeta = calculate_z_vector(theta_high, O2_High_nominal_param, O2_High_cov_matrix, O2_High_max_zeta, O2_High_dm_zeta)

    # Species 5: H2O
    modified_H2O_Low_dm_zeta = calculate_z_vector(theta_low, H2O_Low_nominal_param, H2O_Low_cov_matrix, H2O_Low_max_zeta, H2O_Low_dm_zeta)
    modified_H2O_High_dm_zeta = calculate_z_vector(theta_high, H2O_High_nominal_param, H2O_High_cov_matrix, H2O_High_max_zeta, H2O_High_dm_zeta)

    # Species 6: CO
    modified_CO_Low_dm_zeta = calculate_z_vector(theta_low, CO_Low_nominal_param, CO_Low_cov_matrix, CO_Low_max_zeta, CO_Low_dm_zeta)
    modified_CO_High_dm_zeta = calculate_z_vector(theta_high, CO_High_nominal_param, CO_High_cov_matrix, CO_High_max_zeta, CO_High_dm_zeta)

    # Species 7: C
    modified_C_Low_dm_zeta = calculate_z_vector(theta_low, C_Low_nominal_param, C_Low_cov_matrix, C_Low_max_zeta, C_Low_dm_zeta)
    modified_C_High_dm_zeta = calculate_z_vector(theta_high, C_High_nominal_param, C_High_cov_matrix, C_High_max_zeta, C_High_dm_zeta)

    # Species 8: HCO
    modified_HCO_Low_dm_zeta = calculate_z_vector(theta_low, HCO_Low_nominal_param, HCO_Low_cov_matrix, HCO_Low_max_zeta, HCO_Low_dm_zeta)
    modified_HCO_High_dm_zeta = calculate_z_vector(theta_high, HCO_High_nominal_param, HCO_High_cov_matrix, HCO_High_max_zeta, HCO_High_dm_zeta)

    # Species 9: OH*
    modified_OHstar_Low_dm_zeta = calculate_z_vector(theta_low, OHstar_Low_nominal_param, OHstar_Low_cov_matrix, OHstar_Low_max_zeta, OHstar_Low_dm_zeta)
    modified_OHstar_High_dm_zeta = calculate_z_vector(theta_high, OHstar_High_nominal_param, OHstar_High_cov_matrix, OHstar_High_max_zeta, OHstar_High_dm_zeta)

    # Species 10: H
    modified_H_Low_dm_zeta = calculate_z_vector(theta_low, H_Low_nominal_param, H_Low_cov_matrix, H_Low_max_zeta, H_Low_dm_zeta)
    modified_H_High_dm_zeta = calculate_z_vector(theta_high, H_High_nominal_param, H_High_cov_matrix, H_High_max_zeta, H_High_dm_zeta)

    
    # --- Step 3: Concatenate the 20 modified variables into a single 100-element vector ---
    modified_dm_zeta_vector = np.concatenate([
        modified_AR_Low_dm_zeta, modified_AR_High_dm_zeta,
        modified_H2_Low_dm_zeta, modified_H2_High_dm_zeta,
        modified_O_Low_dm_zeta, modified_O_High_dm_zeta,
        modified_O2_Low_dm_zeta, modified_O2_High_dm_zeta,
        modified_H2O_Low_dm_zeta, modified_H2O_High_dm_zeta,
        modified_CO_Low_dm_zeta, modified_CO_High_dm_zeta,
        modified_C_Low_dm_zeta, modified_C_High_dm_zeta,
        modified_HCO_Low_dm_zeta, modified_HCO_High_dm_zeta,
        modified_OHstar_Low_dm_zeta, modified_OHstar_High_dm_zeta,
        modified_H_Low_dm_zeta, modified_H_High_dm_zeta
    ])    
    # Return the final 100-entry modified vector
    return modified_dm_zeta_vector    

def create_theta_matrix(T_vector):
    """
    Creates the 5x5 temperature basis matrix (Theta) from a vector of 5 temperatures.
    The structure for each row i is [1, T_i, T_i^2, T_i^3, T_i^4].
    
    Args:
        T_vector (np.ndarray): A vector of 5 temperatures (e.g., in Kelvin).

    Returns:
        np.ndarray: The 5x5 Theta matrix.
    """
    T_vector = np.asarray(T_vector)
    if len(T_vector) != 5:
        raise ValueError("T_vector must be a numpy array of length 5.")

    # Stack the powers of T_vector horizontally
    # T**0 is [1, 1, 1, 1, 1]
    # T**1 is [t1, t2, t3, t4, t5]
    # ...
    # np.power(T_vector[:, np.newaxis], np.arange(5)) transposes this for the correct structure
    # where each row is [1, t, t^2, t^3, t^4]
    Theta_matrix = np.row_stack([
        np.ones_like(T_vector),
        T_vector,
        T_vector**2,
        T_vector**3,
        T_vector**4
    ])
    #print(Theta_matrix)
    return Theta_matrix



def calculate_z_vector(
    theta,
    nominal_param,
    cov_matrix,
    max_zeta,
    dm_zeta
):
    """
    Calculates the normalized uncertainty vector 'z' (length 5) by reversing 
    the thermodynamic parameter perturbation process (Formula 30).
    
    The formula is applied element-wise across 5 temperature points:
    
    z_i = (Target_Cp_i - Nominal_Cp_i) / (Max_Cp_i - Nominal_Cp_i)

    Args:
        theta (np.ndarray): Vector of 5 temperatures [t1, t2, t3, t4, t5].
        AR_Low_nominal_param (np.ndarray): The 5 nominal polynomial coefficients (original).
        AR_Low_cov_matrix (np.ndarray): The 5x5 covariance matrix.
        AR_Low_max_zeta (np.ndarray): The 5 maximum perturbation factors (zeta_max).
        AR_Low_dm_zeta (np.ndarray): The 5 target perturbation factors (zeta).

    Returns:
        np.ndarray: The resulting 'z' vector of length 5, where each entry is
                    the normalized deviation factor z_i.
    """
    # --- 1. Construct the Theta Matrix (5x5) ---
    Theta = create_theta_matrix(theta)

    # --- 2. Calculate the Scaled Coefficient Vectors (Length 5) ---
    # These represent the 'perturbed' set of coefficients used in the Cp/R calculation.
    
    # Term 1: Nominal Effect (what the Cp/R would be with original params)
    scaled_nominal = nominal_param

    # Term 2: Maximum Effect (what the Cp/R would be with max positive perturbation)
    scaled_max = np.dot(cov_matrix, max_zeta)

    # Term 3: Target Effect (what the Cp/R is with the target dm_zeta params)
    scaled_target = np.dot(cov_matrix, dm_zeta)

    # --- 3. Calculate the Cp/R values at all 5 points (Length 5 Vectors) ---
    # This is Theta_i^T * (scaled_coefficient_vector) for all 5 rows/points.
    # In numpy, Theta @ vector gives all 5 dot products efficiently.

    # Denominator baseline (Nominal Cp/R at all 5 temperatures)
    Cp_Nominal = np.dot(Theta.T, scaled_nominal)
    
    # Denominator max (Max Cp/R at all 5 temperatures)
    Cp_Max = np.dot(Theta.T, scaled_max)

    # Numerator target (Target Cp/R at all 5 temperatures)
    Cp_Target = np.dot(Theta.T, scaled_target)

    # --- 4. Calculate Numerator and Denominator Vectors ---
    # These are element-wise subtractions across the 5 points.
    
    # Numerator: Actual Deviation from Nominal
    numerator = Cp_Target
    
    # Denominator: Maximum Possible Deviation from Nominal
    denominator = Cp_Max 

    # --- 5. Final Calculation of Z ---
    # Perform element-wise division. Handle potential division by zero.
    z_vector = np.divide(
        numerator, 
        denominator, 
        out=np.zeros_like(numerator), 
        where=denominator!=0
    )
    
    return z_vector    
# =================================================================================
# MAIN EXECUTION
# =================================================================================

def main():
    """
    Main function to read the input CSV, process each row using process_dm_zeta,
    and save the results to a new CSV file.
    """
    # 1. Define the input and output file names
    input_file = "/data2/RANA/FINAL_RUN_7nc/FULL_RUN_1/DesignMatrix.csv"  # <-- CHANGE THIS TO YOUR INPUT FILE NAME
    output_file = "dm_Z_values.csv"
    
 
    print(f"Reading input file: {input_file}...")
    try:
        # Assuming your CSV does not have a header or index
        df_input = pd.read_csv(input_file, header=None)
        rows, cols = df_input.shape
        print(f"The shape of the CSV file is: **{rows}** rows x **{cols}** columns.")
        #raise AssertionError
    except Exception as e:
        print(f"Error reading the input file: {e}")
        return

    # Check if the data size is correct
    if df_input.shape[1] != 100:
        print(f"Error: Input file must have 100 columns (entries per row), but found {df_input.shape[1]}.")
        return

    # 3. Apply the transformation function to every row
    print(f"Processing {df_input.shape[0]} rows...")
    
    # Use the apply method along axis 1 (rows) for efficient row-wise processing.
    # The result will be a pandas Series where each entry is the 100-element modified array.
    df_output_series = df_input.apply(process_dm_zeta, axis=1)

    # 4. Convert the Series of arrays back into a DataFrame
    # This efficiently stacks the 100-element arrays back into a DataFrame.
    df_output = pd.DataFrame(df_output_series.tolist())

    # 5. Save the modified data to the new CSV file
    print(f"Saving modified data to: {output_file}...")
    df_output.to_csv(output_file, header=False, index=False)

    print("\n✅ Processing complete!")
    print(f"Input rows: {df_input.shape[0]}. Output rows: {df_output.shape[0]}.")
    print(f"Output saved to: {output_file}")

if __name__ == "__main__":
    main()

 
