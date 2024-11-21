import cantera as ct
import numpy as np
from pyked import ChemKED
from urllib.request import urlopen
import yaml,time
import pandas as pd
def find_nearest(array, value):
    array = np.asarray(np.abs(np.log(array)))
    value = np.abs(np.log(value))
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

def fast_nearest_interp(xi, x, y):
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]

def ignitionDelay(df,pList, species,cond="max",specified_value ="None;None",exp_conc = "None"):
	if cond == "max":
		tau = df[species].idxmax()
	elif cond == "onset" and species != "p":
		time = df.index.to_numpy()
		conc = (df[species]).to_numpy()
		dtime = np.diff(df.index.to_numpy())
		dconc = np.diff((df[species]).to_numpy())
		slope = dconc/dtime
		intercept = int(np.diff(slope).argmax())
		tau = df.index.to_numpy()[intercept]
	elif cond == "onset" and species == "p":
		time = df.index.to_numpy()
		conc = np.asarray(pList)
		dtime = np.diff(df.index.to_numpy())
		dconc = np.diff(conc)
		slope = dconc/dtime
		intercept = int(np.diff(slope).argmax())
		tau = df.index.to_numpy()[intercept]
	elif cond == "dt-max" and species == "p":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(np.asarray(pList))
		slope = conc/time
		index = int(slope.argmax())
		tau = df.index.to_numpy()[index]
	elif cond == "dt-max" and species != "p":
		time = np.diff(df.index.to_numpy())
		conc = np.diff(df[species].to_numpy())
		slope = conc/time
		index = int(slope.argmax())
		tau = df.index.to_numpy()[index]
	elif cond == "specific":
		if specified_value.split(";")[0] == None:
			raise Assertionerror("Input required for specified_value in ignition delay")
		else:
			target = float(specified_value.split(";")[0])
			unit = specified_value.split(";")[1]
			if unit == "molecule":
				avogadro = 6.02214E+23
				target = target/avogadro
				molecular_wt = gas.atomic_weight(species)
				target = molecular_wt*target	
				index,value = find_nearest(df[species],target)
				tau = df.index.to_numpy()[index[0]]	
			
			elif unit == "molecule/cm3":
				unit_conversion = 10**(6)
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				avogadro = 6.02214E+23
				target = (target/avogadro)*unit_conversion
				if exp_conc != "":
					exp_conc = (exp_conc/avogadro)*unit_conversion
					conc = conc*exp_conc
			
				f = CubicSpline(time,conc, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))			
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
				
			elif unit == "mole/cm3":
				time = df.index.to_numpy()
				conc = df[species].to_numpy()
				if exp_conc != "":
					conc = conc*exp_conc
				
				f = CubicSpline(time,conc, bc_type='natural')
				time_new = np.arange(min(time),max(time),10**(-8))
				conc_new = f(time_new)
				tau = fast_nearest_interp(target,conc_new,time_new)
	return tau

class VolumeProfile(object):
    def __init__(self, keywords):
        self.time = np.array(keywords['vproTime'])
        self.volume = np.array(keywords['vproVol'])/keywords['vproVol'][0]
        self.velocity = np.diff(self.volume)/np.diff(self.time)
        self.velocity = np.append(self.velocity, 0)
    def __call__(self, t):
        if t < self.time[-1]:
            prev_time_point = self.time[self.time <= t][-1]
            index = np.where(self.time == prev_time_point)[0][0]
            return self.velocity[index]
        else:
            return 0

#rcm_link = 'https://raw.githubusercontent.com/pr-omethe-us/PyKED/master/pyked/tests/testfile_rcm.yaml'
#with urlopen(rcm_link) as response:
#    testfile_rcm = yaml.safe_load(response.read())

#ck = ChemKED(dict_input=testfile_rcm,skip_validation=True)
#dp = ck.datapoints[0]

#T_initial = dp.temperature.to('K').magnitude
#P_initial = dp.pressure.to('Pa').magnitude
#X_initial = dp.get_cantera_mole_fraction()

T_initial = 297.4
P_initial = 127722.82894736841
X_initial = {'H2':1.2500e-01, 'O2':6.2500e-02, 'N2':1.8125e-01, 'Ar':6.3125e-01}


gas = ct.Solution('gri30.xml')
gas.TPX = T_initial, P_initial, X_initial

reac = ct.IdealGasReactor(gas)
env = ct.Reservoir(ct.Solution('gri30.xml'))


df = pd.read_csv("volume_profile.csv")
exp_time = df["time(s)"]
exp_volume = df["volume(cm3)"]


stateVariableNames = [reac.component_name(item) for item in range(reac.n_vars)]
timeHistory = pd.DataFrame(columns=stateVariableNames)
#exp_time = dp.volume_history.time.magnitude
#try:
#	exp_volume = dp.volume_history.volume.magnitude
#except:
#	exp_volume = dp.volume_history.quantity.magnitude

#string = "time(s),volume(cm3)\n"
#for i in range(len(exp_time)):
#	string+=f"{exp_time[i]},{exp_volume[i]}\n"
#with open("volume_profile.csv","w") as volume:
#	volume.write(string)

	
keywords = {'vproTime': exp_time, 'vproVol': exp_volume}
ct.Wall(reac, env, velocity=VolumeProfile(keywords));

netw = ct.ReactorNet([reac])
netw.max_time_step = np.min(np.diff(exp_time))

#time = []
#temperature = []
#pressure = []
#volume = []
#mass_fractions = []


t = 0
pressureList = []
counter = 1;
while(gas.T < 1800):
    t = netw.step()
    if (counter%1 == 0):       	
        timeHistory.loc[t] = netw.get_state().astype('float64')
        pressureList.append(gas.P)
    counter+=1

#while netw.time < 0.05:
#    time.append(netw.time)
#    temperature.append(reac.T)
#    pressure.append(reac.thermo.P)
#    volume.append(reac.volume)
#    mass_fractions.append(reac.Y)
#    netw.step()

tic = time.time() 
tau = ignitionDelay(timeHistory,pressureList,"p","onset","None",)
tok = time.time()
tau_file = open("output/tau.out",'w')
tau_file.write("#T(K)    tau(us)"+"\n"+"{}    {}".format(T_initial,tau*(10**6)))
tau_file.close()   
 
#timeHistory.to_csv("time_history.csv")
#For plotting purpose

#import matplotlib.pyplot as plt
#plt.figure()
#plt.plot(time, pressure)
#plt.ylabel('Pressure [Pa]')
#plt.xlabel('Time [s]');
#plt.show()

#plt.figure()

#plt.plot(exp_time, exp_volume/exp_volume[0], label='Experimental volume', linestyle='--')
#plt.plot(time, volume, label='Simulated volume')
#plt.legend(loc='best')
#plt.ylabel('Volume [m^3]')
#plt.xlabel('Time [s]');
#plt.show()

