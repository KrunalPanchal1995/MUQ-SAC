import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Reading simulation data
df_sim = pd.read_csv("p4_.csv")
df_sim = pd.DataFrame(df_sim)
dT = 1000/df_sim["T"].to_numpy()
Sim = df_sim["Tig(us)"].to_numpy()/1000

#Reading observed data
df_obs = pd.read_csv("P_4_exp.csv")
df_obs = pd.DataFrame(df_obs)
dT_obs = df_obs["1000/T"].to_numpy()
Obs = df_obs["Obs(ms)"].to_numpy()

#Reading simulation data from the literature Lele at al 2018
df_sim_lele = pd.read_csv("P4_sim_lele.csv")
df_sim_lele = pd.DataFrame(df_sim_lele)
dT_sim_lele = df_sim_lele["1000/T"].to_numpy()
Sim_lele = df_sim_lele["Sim(ms)"].to_numpy()
#For plotting the data

fig = plt.figure()
plt.plot(dT,Sim,"-",label="MB ver2.0 (cantera)")
plt.plot(dT_obs,Obs,"ro",label="Lele et al. 2018")
plt.plot(dT_sim_lele,Sim_lele,"r--",label="Lele simulation")
plt.legend()
plt.yscale("log")
plt.xlabel("1000/T")
plt.ylabel("Ignition delay time (ms)")
plt.title(r"P = 4 bar, $\phi$ = 0.25")
#plt.show()
plt.savefig("P_4.pdf")
