import numpy as np
import pandas as pd

file_ = "/home/krithika/MUQ-SAC/Database/T124MBZ/TMB/rcm_profiles/RCM_TBM_100_PHI_1_P_30/T_697.csv"
df = pd.read_csv(file_, delim_whitespace=True, header=None)
df.columns = ["Tag", "Time", "Volume"]
print(df["Time"].to_numpy())
