import xml.etree.ElementTree as ET
import scipy as sp
import numpy as np
import yaml
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import os,re
import matplotlib.pyplot as plt
from scipy.linalg import block_diag
from pyDOE2 import *

class UQ(object):
    def __init__(self,data):
        self.data = data
        self.tempData = data["temperatures"]
        self.unsrtData = data["uncertainties"]
        self.temperatures = []
        self.uncertainties = []
        self.foundCov = None
        self.getArrays(flag=data["flag"])
        self.ArrheniusParams = np.array([float(data["Arrhenius"]["A"]),float(data["Arrhenius"]["n"]),float(data["Arrhenius"]["Ea"])])
        if "log" in self.data["Arrhenius"]["flag"]:
            pass
        else:
            self.ArrheniusParams = np.array([np.log(self.ArrheniusParams),self.ArrheniusParams[1],self.ArrheniusParams[2]])
        if "Covariance" in data:
            self.test = []
            self.text_cov = data["Covariance"]
            for i in range(3):
                self.test.append(np.asfarray(self.text_cov["mat"+str(i+1)].split(","),float))
            #print(np.matrix(self.test))
            self.foundCov = True
        if "test_samples" in data:
            file_name = data["test_samples"]
            sampling_file = open(file_name,"r").readlines()
            self.Nagy_arrhenius_samples = []
            for sample in sampling_file[1:]:
                self.Nagy_arrhenius_samples.append(np.asfarray(sample.strip("''").strip("\n").split()[:3],float))
        self.Theta = np.array([self.temperatures/self.temperatures,np.log(self.temperatures),-1/(self.temperatures)])
        self.M = 3.0/np.log(10.0)
        self.guess = np.array([-10.0,-10.0,0.5,200.0,10.0,5.0,1,1])
        self.guess_z = np.array([1,1,1])
        
    def getArrays(self,flag=False):
        if flag==True:
            nBounds = int(self.tempData["nB"])
            try:
                nB_unsrt = int(self.unsrtData["nB"])
                if nB_unsrt != nBounds:
                    raise AssertionError("Give proper number of bounds")
            except:
                print("Uncertainty check completer !!")
            for i in range(nBounds):
                temp = []
                unsrt = []
                a = self.tempData["b"+str(i+1)]["a"]
                b = self.tempData["b"+str(i+1)]["b"]
                temp = list(np.linspace(a,b,100))
                unsrt = list(self.unsrtData["u"+str(i+1)]*np.ones(len(temp)))
                self.temperatures.extend(temp)
                self.uncertainties.extend(unsrt)
                
            self.temperatures = np.asarray(self.temperatures)
            self.uncertainties = np.asarray(self.uncertainties)
        else:
            unsrt_file = open(self.data["unsrtFile"],"r").readlines()
            unsrtData = [np.asfarray(i.strip("\n").strip("''").split(","),float) for i in unsrt_file]
            self.temperatures = np.asarray([i[0] for i in unsrtData])
            self.uncertainties = np.asarray([i[1] for i in unsrtData])

    def getUncertFunc(self,L):
        #func = [self.M*np.linalg.norm(np.dot(L.T,i)) for i in self.Theta.T]
        func = self.M*np.sqrt((L[0][0]*self.Theta[0]+L[1][0]*self.Theta[1]+L[2][0]*self.Theta[2])**2+(L[1][1]*self.Theta[1]+L[2][1]*self.Theta[2])**2+(L[2][2]*self.Theta[2])**2)
        return np.asarray(func)
    
    def getZetaUnsrtFunc(self,L,z):
        #func = [self.M*(i.T.dot(L.dot(z))) for i in self.Theta.T]
        #print(func)
        func = self.M*((L[0][0]*z[0])*self.Theta[0]+(L[1][0]*z[0]+L[1][1]*z[1])*self.Theta[1]+(L[2][0]*z[0]+L[2][1]*z[1]+L[2][2]*z[2])*self.Theta[2])
        #print(func)
        return np.asarray(func)
    
    def obj_func_eigen(self,guess):
        z = guess
        cov = np.array([[z[0],z[1],z[2]],[z[3],z[4],z[5]],[z[6],z[7],z[8]]]);#cholesky lower triangular matric
        #self.unsrtFunc = [self.M*np.linalg.norm(np.dot(cov.T,i.T)) for i in self.Theta.T]        
        f = (self.uncertainties - self.M*self.getUncertFunc(cov))/(self.uncertainties/self.M)
        obj = np.dot(f,f)
        return obj
    
    def obj_func(self,guess):
        z = guess
        cov = np.array([[z[0],0,0],[z[1],z[2],0],[z[3],z[4],z[5]]]);#cholesky lower triangular matric
        #f = np.zeros(len(self.temperatures))
        #self.unsrtFunc = [self.M*(np.dot(i.T,np.dot(cov,np.dot(cov.T,i)))) for i in self.Theta.T]
        #self.unsrtFunc = [self.M*np.linalg.norm(np.dot(cov.T,i.T)) for i in self.Theta.T]        
        f = (self.uncertainties - self.M*self.getUncertFunc(cov))/(self.uncertainties/self.M)
        obj = np.dot(f,f)
        return obj
    
    def const_func(self,guess):
        self.z = guess
        f = np.zeros(len(self.temperatures))
        f = self.uncertainties - self.getUncertFunc()
        return np.amin(f)
    
    def obj_func_zeta(self,guess):
        M = self.M
        T = self.temperatures
        cov = self.L
        Theta = np.array([T/T,np.log(T),-1/(T)])
        f = (self.uncertainties-self.getZetaUnsrtFunc(cov,guess))
        obj = np.dot(f,f)
        return obj
    
    
    def const_func_zeta_1(self,z):
        M = self.M
        T = self.temperatures
        cov = self.L
        Theta = np.array([T/T,np.log(T),-1/(T)])
        normTheta = np.sqrt((T/T)**2 + (np.log(T))**2 + (1/T)**2)
        unsrtFunc = M*np.linalg.norm(np.dot(cov.T,Theta))
        uncorrFunc = np.linalg.norm(np.dot(cov,guess))
        QLT = np.asarray(np.dot(Theta.T,np.dot(cov,guess))).flatten()
        f = (self.uncertainties - self.getZetaUnsrtFunc(cov,z))
        #print(np.dot(np.transpose(theta),np.dot(L1,eta)))
        obj = np.dot(f,f)
        return np.amin(f)

    def const_func_zeta_2(self,z):
        M = self.M
        T = self.temperatures[1:-1]
        cov = self.L
        Theta = np.array([T/T,np.log(T),-1/(T)])
        QtLZ = np.asarray([self.M*(i.dot(cov.dot(z))) for i in Theta.T])
        f = (self.uncertainties[1:-1]-QtLZ)
        return np.amin(f)
    
    def const_1_typeC_Zeta(self,z):
        M = self.M
        #T = self.temperatures
        T = self.temperatures[0]
        cov = self.L
        Pmin = self.P_min
        P = self.ArrheniusParams
        Theta = np.array([T/T,np.log(T),-1/(T)])
        k_min = Theta.dot(Pmin)
        QtLZ = self.M*(Theta.T.dot(cov.dot(z)))
        f = Theta.dot(P)-QtLZ
        return k_min - f
    
    def const_2_typeC_Zeta(self,z):
        M = self.M
        #T = self.temperatures
        T = self.const_T
        cov = self.L
        Pmax = self.P_max
        P = self.ArrheniusParams
        Theta = np.array([T/T,np.log(T),-1/(T)])
        k_max = Theta.dot(Pmax)
        QtLZ = self.M*(Theta.T.dot(cov.dot(z)))
        f = (Theta.dot(P)-QtLZ)
        return k_max - f
    
    def const_2_typeA_Zeta(self,z):
        M = self.M
        T = self.temperatures[-1]
        cov = self.L
        Pmax = self.P_max
        P = self.ArrheniusParams
        Theta = np.array([T/T,np.log(T),-1/(T)])
        k_max = Theta.dot(Pmax)
        QtLZ = self.M*(Theta.T.dot(cov.dot(z)))
        f = (Theta.dot(P)-QtLZ)
        return k_max - f
    
    def const_2_typeA_Zeta(self,z):
        M = self.M
        T = self.temperatures[-1]
        cov = self.L
        Pmax = self.P_max
        P = self.ArrheniusParams
        Theta = np.array([T/T,np.log(T),-1/(T)])
        k_max = Theta.dot(Pmax)
        QtLZ = self.M*(Theta.T.dot(cov.dot(z)))      
        f = (Theta.dot(P)-QtLZ)
        return k_max - f
    
    def const_3_typeC_Zeta(self,z):
        M = self.M
        #T = self.temperatures
        T = self.temperatures[-1]
        cov = self.L
        Pmin = self.P_min
        P = self.ArrheniusParams
        Theta = np.array([T/T,np.log(T),-1/(T)])
        k_min = Theta.dot(Pmin)
        QtLZ = self.M*(Theta.T.dot(cov.dot(z)))  
        f = (Theta.dot(P)-QtLZ)
        return k_min - f
    
    
    def getCovariance(self,flag = False):
        if flag == True:
            constraints = {'type': 'ineq', 'fun': self.const_func }
            self.const = [constraints]
            self.solution = minimize(self.obj_func,self.guess,constraints=self.cons)
        else:
            self.solution = minimize(self.obj_func,self.guess,method="Nelder-Mead")
        self.L = np.array([[self.solution.x[0],0,0],[self.solution.x[1],self.solution.x[2],0],[self.solution.x[3],self.solution.x[4],self.solution.x[5]]]);#cholesky lower triangular matric
        cov1 = self.L
        cov2 = np.dot(self.L,self.L.T)
        #print(cov2)
        #print(np.exp(self.ArrheniusParams[0]),self.ArrheniusParams[1],self.ArrheniusParams[2])
        #print(self.temperatures[0],self.temperatures[-1])
        D,Q = np.linalg.eigh(cov2)
        self.A = Q.dot(sp.linalg.sqrtm(np.diag(D)))
        return self.L
    
    def getUnCorrelated(self,flag = False):
        if flag == True:
            #con1 = {'type': 'ineq', 'fun': self.const_func_zeta_1}
            #con2 = {'type': 'ineq', 'fun': self.const_func_zeta_2}
            con3 = {'type': 'eq', 'fun': self.const_nmax}
            con4 = {'type': 'eq', 'fun': self.const_nmin}
            self.const_zeta = [con3,con4]
            zeta = minimize(self.obj_func_zeta,self.guess_z,constraints=self.const_zeta)
        else:
            zeta = minimize(self.obj_func_zeta,self.guess_z,method="Nelder-Mead")
        return zeta
    
    def getConstrainedUnsrtZeta(self,flag=False):
        if flag == True:
            con1 = {'type': 'eq', 'fun': self.const_1_typeC_Zeta}
            con2 = {'type': 'eq', 'fun': self.const_2_typeC_Zeta}
            con3 = {'type': 'eq', 'fun': self.const_3_typeC_Zeta}
            self.const_zeta = [con1,con2,con3]
            zeta_list = []
            obj_val = []
            alpha = []
            n = []
            epsilon = []
            
            for i,T in enumerate(self.temperatures):
                if T<self.temperatures[-1]-300 and T> self.temperatures[0]+300:
                    self.const_T = T
                    zeta = minimize(self.obj_func_zeta,self.guess_z,method="SLSQP",constraints=self.const_zeta)
                    #zeta = minimize(self.obj_curved_zeta,self.guess_z,method="Nelder-Mead")
                    zeta_list.append(zeta.x)
                    alpha.append(zeta.x[0])
                    n.append(zeta.x[1])
                    epsilon.append(zeta.x[2])
                    obj_val.append(abs(self.obj_func_zeta(zeta.x))+abs(self.const_1_typeC_Zeta(zeta.x))+abs(self.const_2_typeC_Zeta(zeta.x))+abs(self.const_3_typeC_Zeta(zeta.x)))
            alpha_square = [i**2 for i in alpha]
            n_square = [i**2 for i in n]
            epsilon_square = [i**2 for i in epsilon]
            index = obj_val.index(min(obj_val))
            
            return [zeta_list[index]],index,np.array([alpha[alpha_square.index(max(alpha_square))],n[n_square.index(max(n_square))],epsilon[epsilon_square.index(max(epsilon_square))]])
        else:
            con1 = {'type': 'eq', 'fun': self.const_1_typeC_Zeta}
            con2 = {'type': 'eq', 'fun': self.const_2_typeA_Zeta}
            con3 = {'type': 'ineq', 'fun': self.const_func_zeta_2}
            self.const_zeta = [con1,con2,con3]
            zeta = minimize(self.obj_func_zeta,self.guess_z,constraints=self.const_zeta)
            return zeta
        #return zeta_list
    
    def getUnsrtData(self,tag):
        T = self.temperatures
        P = self.ArrheniusParams
        self.cov = self.getCovariance()
        self.zeta = self.getUnCorrelated(flag=False)
        self.unsrtFunc = self.getUncertFunc(self.cov)
        self.zetaUnsrt = self.getZetaUnsrtFunc(self.cov,self.zeta.x)
       
        self.P_max = P + np.asarray(np.dot(self.cov,self.M*zeta)).flatten();
        self.P_min = P - np.asarray(np.dot(self.cov,self.M*zeta)).flatten();
        self.kmax = self.Theta.T.dot(self.P_max)
        self.kmin = self.Theta.T.dot(self.P_min)
        self.kappa = self.Theta.T.dot(P)
        self.zeta_curved_type_A,index,zeta_lim = self.getConstrainedUnsrtZeta(flag=True)
        self.zeta_curved_type_C = self.getConstrainedUnsrtZeta(flag=False)
        """
        fig, axs = plt.subplots(2, 1, figsize=(15,20))

        for zeta in self.zeta_curved_type_A:
            temp_unsrtFunc = np.asarray([self.M*np.dot(i.T,np.dot(self.cov,zeta)) for i in self.Theta.T])
            axs[0].plot(self.temperatures,temp_unsrtFunc,'k--')
            axs[0].plot(self.temperatures,-temp_unsrtFunc,'k--')
            
        axs[0].set_xlabel('Temperature (K)')
        axs[0].set_ylabel('Uncertainity ($f$)')
        temp_unsrtFunc = np.asarray([np.dot(i.T,np.dot(self.cov,self.M*self.zeta.x)) for i in self.Theta.T])
        temp_unsrtFunc_C = np.asarray([np.dot(i.T,np.dot(self.cov,self.M*self.zeta_curved_type_C.x)) for i in self.Theta.T])

        axs[0].plot(self.temperatures,temp_unsrtFunc_C,'y--')
        axs[0].plot(self.temperatures,-temp_unsrtFunc_C,'y--')
        if self.foundCov == True:
            unknownZeta = []
            alpha = []
            n = []
            epsilon = []
            self.testcov = np.matrix(self.test)
            L = np.linalg.cholesky(self.testcov)
            self.test_unsrtFunc = [self.M*np.linalg.norm(np.dot(L.T,i)) for i in self.Theta.T]
        #axs[0].plot(self.temperatures,self.unsrtFunc,'r--',label='present study (MUQ)')
        #axs[0].plot(self.temperatures,-self.unsrtFunc,'r--')
        axs[0].plot(self.temperatures,self.zetaUnsrt,'b-',label='present study (zeta) (MUQ)')
        axs[0].plot(self.temperatures,-self.zetaUnsrt,'b-')
        axs[0].set_ylim(-2*max(self.unsrtFunc),2*max(self.unsrtFunc))
        axs[0].plot(self.temperatures,self.uncertainties,'go',label='Exp. data')
        
        
        ###################
               
        
        #fig = plt.figure()
        axs[1].set_xlabel(r"$1/T(1/K)$")
        axs[1].set_ylabel(r"Rate constant ($\{\kappa \}= \log(k)$)}$")
        axs[1].plot((1/T),self.kappa,"-",color="#15B01A",linewidth=1.5,label='$\kappa_o+f_{prior}(T)$')
        axs[1].plot((1/T),self.kmax,":",color ='#030764',linewidth=1.7,label=r'$\kappa_{max}\,(\phi_{\zeta} = f(T)\minus \theta^T L\vec{\zeta}$)')
        axs[1].plot((1/T),self.kmin,':',color='#580F41',linewidth=1.7,label=r'$\kappa_{min}\,(\phi_{\zeta} = f(T)\minus \theta^T L\vec{\zeta}$)')
        count = 0
        beta_selected = []
        """
        zeta = np.array([[self.zeta.x[0],self.zeta.x[1],self.zeta.x[2]],[-self.zeta_curved_type_A[0][0],-self.zeta_curved_type_A[0][1],-self.zeta_curved_type_A[0][2]],[self.zeta_curved_type_C.x[0],self.zeta_curved_type_C.x[1],self.zeta_curved_type_C.x[2]]])
        self.zeta_matrix = np.matrix(zeta)
        """
        beta_list = bbdesign(3)
        for beta in beta_list:
            samples = P + np.asarray(np.dot(self.cov,zeta.T.dot(beta.T).T)).flatten();
            k = self.Theta.T.dot(samples)
            diff = []
            diff_max = []
            for i,ele in enumerate(self.kappa):
                diff.append(np.abs(k[i]-self.kappa[i]))
                diff_max.append(np.abs(self.kappa[i]-self.kmin[i]))
           
            if max(diff)>max(diff_max):
                axs[1].plot((1/T),k,"k-",linewidth=0.35)
                continue
            else:
                axs[1].plot((1/T),k,"r-",linewidth=0.35)
                count +=1
                beta_selected.append(beta)
        for i in range(10):
            var = np.array([0.33,0.33,0.33])
            mean = np.zeros(3)
            beta = np.asarray(np.random.normal(mean,var,3))
            samples = P + np.asarray(np.dot(self.cov,zeta.T.dot(beta.T).T)).flatten();
            k = self.Theta.T.dot(samples)
            diff = []
            diff_max = []
            for i,ele in enumerate(self.kappa):
                diff.append(np.abs(k[i]-self.kappa[i]))
                diff_max.append(np.abs(self.kappa[i]-self.kmin[i]))
           
            if max(diff)>max(diff_max):
                #axs[1].plot((1/T),k,"k-",linewidth=0.15)
                continue
            else:
                #axs[1].plot((1/T),k,"r-",linewidth=0.15)
                count +=1
                beta_selected.append(beta)
        
        
        if self.foundCov == True:
            for i in range(20):
                samples = self.Nagy_arrhenius_samples[i]
                
                k = self.Theta.T.dot(samples)
                #axs[1, 1].plot((1/T),k,"k-",linewidth=0.35)
            
        for ax in axs.flat:
            ax.label_outer()
        plt.legend(fancybox=True, framealpha=0)
        plt.savefig("Sampling_plots/rate_constant_sampling"+str(tag)+".svg",format='svg', dpi=1200,bbox_inches="tight")
        plt.show();
    """
   

# In[23]:


#import yaml file
#Input_file = open("Input.yaml",'r').read()
#data = yaml.safe_load(Input_file)







