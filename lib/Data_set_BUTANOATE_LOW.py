import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import mpl_toolkits.axisartist as AA
from scipy.interpolate import make_interp_spline

plt.style.use('seaborn-white')
import pandas as pd
import os
import csv

cm = 1/2.54 


files = [i for i in os.listdir() if i.endswith("csv")]

dataset_Tig = []
dataset_Fls_phi = []
dataset_Fls_P = []
plot_data_Fls_phi = {}
plot_data_Fls_P = {}
plot_data_Tig = {}
datasets = []


for file_ in files:
	
	if "x1" or "ign" in file_:
		#print(file_)
		temp = {}
		r = file_.strip(".csv")
		dataset_Tig.append(file_.strip(".csv"))
		data = pd.read_csv(file_,sep=',',engine="python")
		df = pd.DataFrame(data,columns=["1000/T\K^-1","obs","std","Nominal","Opt","PRS"])
		temp["Temp"] = df["1000/T\K^-1"].tolist()
		temp["obs"] = df["obs"].tolist()
		temp["std"] = df["std"].tolist()
		temp["Nominal"] = df["Nominal"].tolist()
		temp["Opt"] = df["Opt"].tolist()
		temp["PRS"] = df["PRS"].tolist()
		plot_data_Tig[r] = temp
		
	if "x2" in file_:	
		#print(file_)
		temp = {}
		r = file_.strip(".csv")
		dataset_Tig.append(file_.strip(".csv"))
		data = pd.read_csv(file_,sep=',',engine="python")
		if "Phi" in data:
			#print(f"Phi,{file_}\n")
			dataset_Fls_phi.append(file_.strip(".csv"))
			df = pd.DataFrame(data,columns=["Phi","obs","std","Nominal","Opt","PRS"])
			temp["Phi"] = df["Phi"].tolist()
			temp["obs"] = df["obs"].tolist()
			temp["std"] = df["std"].tolist()
			temp["Nominal"] = df["Nominal"].tolist()
			temp["Opt"] = df["Opt"].tolist()
			temp["PRS"] = df["PRS"].tolist()
			plot_data_Fls_phi[r] = temp
			
		elif "P[atm]" in data:
			#print(f"P[atm],{file_}\n")
			dataset_Fls_phi.append(file_.strip(".csv"))
			df = pd.DataFrame(data,columns=["P[atm]","obs","std","Nominal","Opt","PRS"])
			temp["Pressure"] = df["P[atm]"].tolist()
			temp["obs"] = df["obs"].tolist()
			temp["std"] = df["std"].tolist()
			temp["Nominal"] = df["Nominal"].tolist()
			temp["Opt"] = df["Opt"].tolist()
			temp["PRS"] = df["PRS"].tolist()
			plot_data_Fls_P[r] = temp
	datasets.append(r)
	
#print(datasets)

datasets = ['ing_Butanoate','ing_Ethyl_Butanoate']


#datasets = ['x10001023', 'x10001053', 'x10001029', 'x10001051','x10001054', 'x10001030', 'x20001004','x20001005','x20001034','x20001035','x20001040','x20001176','x20001177','x20001178']

#datasets = ['x20001140', 'x20001035', 'x20001003', 'x10001062', 'x10001055', 'x10001083', 'x20001182', 'x20001068', 'x20001178', 'x10001023', 'x20001176', 'x20001034', 'x10001053', 'x20001180', 'x10001031', 'x10001029', 'x20001040', 'x10001022', 'x20001181', 'x20001177', 'x10001071', 'x20001179', 'x10001051', 'x10001068', 'x10001049', 'x20001044', 'x10001030', 'x10001080', 'x10001065', 'x10001028']

#print(plot_data_Tig)

"""
Ignition delay plots
"""
T = {}
P = {}
Obs = {}
std_vnt = {}
prior = {}
opt = {}
PRS = {}

for i in plot_data_Tig:
	T[i] = plot_data_Tig[i]["Temp"]
	Obs[i] = plot_data_Tig[i]["obs"]
	std_vnt[i] = plot_data_Tig[i]["std"]
	prior[i] = plot_data_Tig[i]["Nominal"]
	opt[i] = plot_data_Tig[i]["Opt"]
	PRS[i] = plot_data_Tig[i]["PRS"]

for i in plot_data_Fls_P:
	P[i] = plot_data_Fls_P[i]["Pressure"]
	Obs[i] = plot_data_Fls_P[i]["obs"]
	std_vnt[i] = plot_data_Fls_P[i]["std"]
	prior[i] = plot_data_Fls_P[i]["Nominal"]
	opt[i] = plot_data_Fls_P[i]["Opt"]
	
	
for i in plot_data_Fls_phi:
	P[i] = plot_data_Fls_phi[i]["Phi"]
	Obs[i] = plot_data_Fls_phi[i]["obs"]
	std_vnt[i] = plot_data_Fls_phi[i]["std"]
	prior[i] = plot_data_Fls_phi[i]["Nominal"]
	opt[i] = plot_data_Fls_phi[i]["Opt"]
"""
Group_0: ign_phi1_0_p13_5
"""	
diff_23_P = (np.asarray(prior["ing_Butanoate"])-np.asarray(Obs["ing_Butanoate"]))
diff_23_O = (np.asarray(opt["ing_Butanoate"])-np.asarray(Obs["ing_Butanoate"]))#/
diff_23_PRS=(np.asarray(PRS["ing_Butanoate"])-np.asarray(Obs["ing_Butanoate"]))

#np.asarray(Obs["ign_phi1_0_p13_5"])

RMS_23_P = np.sqrt(np.mean(diff_23_P**2))
RMS_23_O = np.sqrt(np.mean(diff_23_O**2))
improvement_23 = (RMS_23_P-RMS_23_O)*100/RMS_23_P
fig,ax = plt.subplots()


from scipy.interpolate import make_interp_spline

#x = np.asarray(T["ing_Butanoate"]).flatten()
#y = np.asarray(prior["ing_Butanoate"]).flatten()

#X_Y_Spline = make_interp_spline(x,y)

# Returns evenly spaced numbers
# over a specified interval.
#X_ = np.linspace(x.min(), x.max(), 500)
#Y_ = X_Y_Spline(X_)
#ax.text(0.8,6000,r'MB-C5H10O2/O$_2$/AR, $\phi$ = 1',ha="left",va="top", fontsize=10, color='black')
#ax.text(0.42,3700,r'Dilution 1:5, N-C7H16$_2$',ha="left",va="top", fontsize=10, color='black')
plt.yscale('log')
#plt.xlabel(r"$1000/T (K^{-1})$")
#plt.ylabel(r"Ignition delay, $\tau$ ($\mu s$)")
#plt.plot(X_,Y_,"r-.",linewidth=1.1,label="Prior simulation")
plt.plot(T["ing_Butanoate"],opt["ing_Butanoate"],"b-",linewidth=1.1,label="Optimized mechanism")

plt.errorbar(T["ing_Butanoate"],(np.asarray(Obs["ing_Butanoate"])),yerr = (np.asarray(std_vnt["ing_Butanoate"])),fmt='ks',ecolor="black",markerfacecolor='black',markeredgecolor='black',markersize=4,capsize=2,elinewidth=0.7,markeredgewidth=0.5,label = f"(a) $P$ = 13.5 atm")#, Prior RMS = {RMS_23_P:.2f}, Optimized RMS = {RMS_23_O:.2f}")

plt.legend(loc="best")
#plt.ylim(50,40000)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', width=1.2,direction = 'in')
ax.tick_params(which='major', length=5)
ax.tick_params(which='minor', length=3, color='k')
plt.savefig("group0.jpeg",bbox_inches="tight")
plt.show()
#print(T)



diff_23_P = (np.asarray(prior["ing_Butanoate"])-np.asarray(Obs["ing_Butanoate"]))
diff_23_O = (np.asarray(opt["ing_Butanoate"])-np.asarray(Obs["ing_Butanoate"]))#/np.asarray(Obs["ign_phi1_0_p13_5"])
diff_23_PRS=(np.asarray(PRS["ing_Butanoate"])-np.asarray(Obs["ing_Butanoate"]))
RMS_23_P = np.sqrt(np.mean(diff_23_P**2))
RMS_23_O = np.sqrt(np.mean(diff_23_O**2))
RMS_23_PRS = np.sqrt(np.mean(diff_23_PRS**2))


diff_29_P = (np.asarray(prior["ing_Ethyl_Butanoate"])-np.asarray(Obs["ing_Ethyl_Butanoate"]))#/np.asarray(Obs["ign_phi2_0_p13_5"])
diff_29_O = (np.asarray(opt["ing_Ethyl_Butanoate"])-np.asarray(Obs["ing_Ethyl_Butanoate"]))#/np.asarray(Obs["ign_phi2_0_p13_5"])
diff_29_PRS = (np.asarray(PRS["ing_Ethyl_Butanoate"])-np.asarray(Obs["ing_Ethyl_Butanoate"]))

RMS_29_P = np.sqrt(np.mean(diff_29_P**2))
RMS_29_O = np.sqrt(np.mean(diff_29_O**2))
RMS_29_PRS = np.sqrt(np.mean(diff_29_PRS**2))
improvement_23 = (RMS_23_P-RMS_23_O)*100/RMS_23_P
improvement_29 = (RMS_29_P-RMS_29_O)*100/RMS_29_P

n_model_prior_23 = np.log(np.asarray(prior["ing_Butanoate"])*10)
n_model_opt_23 = np.log(np.asarray(opt["ing_Butanoate"])*10)
n_model_PRS_23 = np.log(np.asarray(PRS["ing_Butanoate"])*10)
n_exp_23 = np.log(np.asarray(Obs["ing_Butanoate"])*10)
s_exp_23 = np.asarray(std_vnt["ing_Butanoate"])/np.asarray(Obs["ing_Butanoate"])

obj_prior_23 = (n_model_prior_23-n_exp_23)/s_exp_23
o_23_p = 0
for dif in obj_prior_23:
	o_23_p+=dif**2
obj_opt_23 = (n_model_opt_23-n_exp_23)/s_exp_23
o_23_o = 0
for dif in obj_opt_23:
	o_23_o+=dif**2

obj_PRS_23 = (n_model_PRS_23-n_exp_23)/s_exp_23
o_23_PRS =0
for dif in obj_PRS_23:
	o_23_PRS+=dif**2



n_model_prior_29 = np.log(np.asarray(prior["ing_Ethyl_Butanoate"])*10)
n_model_opt_29 = np.log(np.asarray(opt["ing_Ethyl_Butanoate"])*10)
n_model_PRS_29 = np.log(np.asarray(PRS["ing_Ethyl_Butanoate"])*10)
n_exp_29 = np.log(np.asarray(Obs["ing_Ethyl_Butanoate"])*10)
s_exp_29 = np.asarray(std_vnt["ing_Ethyl_Butanoate"])/np.asarray(Obs["ing_Ethyl_Butanoate"])

obj_prior_29 = (n_model_prior_29-n_exp_29)/s_exp_29
o_29_p = 0
for dif in obj_prior_29:
	o_29_p+=dif**2
obj_opt_29 = (n_model_opt_29-n_exp_29)/s_exp_29
o_29_o = 0
for dif in obj_opt_29:
	o_29_o+=dif**2
obj_PRS_29 = (n_model_PRS_29-n_exp_29)/s_exp_29
o_29_PRS = 0
for dif in obj_PRS_29:
	o_29_PRS+=dif**2

fig,ax = plt.subplots()
ax.text(0.8,6000,r'MB-C5H10O2/O$_2$/AR',ha="left",va="top", fontsize=10, color='black')
#ax.text(0.42,3700,r'Dilution 1:5, N-C7H16$_2$',ha="left",va="top", fontsize=10, color='black')
plt.yscale('log')
plt.xlabel(r"$1000/T (K^{-1})$")
plt.ylabel(r"Ignition delay, $\tau$ ($\mu s$)")
plt.plot(T["ing_Butanoate"],prior["ing_Butanoate"],"r-.",linewidth=1.1,label="Prior simulation")
plt.plot(T["ing_Butanoate"],opt["ing_Butanoate"],"b-",linewidth=1.1,label="Optimized mechanism")
plt.plot(T["ing_Butanoate"],PRS["ing_Butanoate"],"b--",linewidth=1.1,label="PRS")

plt.errorbar(T["ing_Butanoate"],(np.asarray(Obs["ing_Butanoate"])),yerr = (np.asarray(std_vnt["ing_Butanoate"])),fmt='ks',ecolor="black",markerfacecolor='black',markeredgecolor='black',markersize=4,capsize=2,elinewidth=0.7,markeredgewidth=0.5,label = f"(a)$\phi$ = 1, $P$ = 30 bar, Prior obj = {o_23_p:.2f}, Optimized obj = {o_23_o:.2f}, PRS = {o_23_PRS:.2f}")
plt.plot(T["ing_Ethyl_Butanoate"],prior["ing_Ethyl_Butanoate"],"r-.",linewidth=1.1)
plt.plot(T["ing_Ethyl_Butanoate"],opt["ing_Ethyl_Butanoate"],"b-",linewidth=1.1)
plt.plot(T["ing_Ethyl_Butanoate"],PRS["ing_Ethyl_Butanoate"],"b--",linewidth=1.1)
plt.errorbar(T["ing_Ethyl_Butanoate"],(np.asarray(Obs["ing_Ethyl_Butanoate"])),yerr =(np.asarray(std_vnt["ing_Ethyl_Butanoate"])),fmt='k.',ecolor="black",markerfacecolor='black',markeredgecolor='black',markersize=8,capsize=2,elinewidth=0.7,
    markeredgewidth=0.5,label = f"(b)$\phi$ = 1, $P$ = 30 bar, Prior obj = {o_29_p:.2f}, Optimized obj = {o_29_o:.2f}, PRS = {o_29_PRS:.2f}")
plt.legend(loc="best")
#plt.ylim(50,40000)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', width=1.2,direction = 'in')
ax.tick_params(which='major', length=5)
ax.tick_params(which='minor', length=3, color='k')
plt.savefig("group1.pdf",bbox_inches="tight")
plt.show()

