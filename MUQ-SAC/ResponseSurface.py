import os
import math
import json
import statistics
import numpy as np
import scipy as sp
from numpy import linalg as LA
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor as MLP
from sklearn.preprocessing import PolynomialFeatures
from scipy.interpolate import InterpolatedUnivariateSpline
from sklearn.datasets import make_friedman2
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.linear_model import HuberRegressor, LinearRegression
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import QuantileRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel

class ResponseSurface(object):
	def __init__(self,xdata,ydata,case,case_index,responseOrder=2):
		self.X = xdata
		self.Y = ydata
		self.case = case
		self.case_index = case_index
		self.order = int(responseOrder)
	
	def create_Neural_Network(self):
		self.regr = MLP(hidden_layer_sizes=(31,30,15),max_iter=1500).fit(self.X,self.Y)
		self.regr_fit = self.regr.predict(self.X)
		#self.regr.score(self.X,self.Y)
		print(self.regr.score(self.X,self.Y))
	def create_SVR_response_surface(self):
		#print(self.Y)
		self.svr = SVR().fit(self.X,self.Y)
		self.svr_yfit = self.svr.predict(self.X)

	def test(self,xTest,yTest):
		self.y_Test_Predict = []
		self.y_Test_simulation = yTest
		for sample in xTest:
			self.y_Test_Predict.append(self.evaluate(sample))
		
		self.error_testing = []
		self.error_testing_relative = []
		for index,sample in enumerate(self.y_Test_simulation):
			self.error_testing.append(float(np.asarray(self.y_Test_Predict)[index])-sample)
			self.error_testing_relative.append((self.y_Test_Predict[index]-sample)/(sample)*100)
		
		self.ytestMaxError = max(self.error_testing_relative)
		self.yTestMeanError = statistics.mean(self.error_testing_relative)
		if self.ytestMaxError > 2:#if maximum error in PRS testing is greater than 2 percentage, the response surface is rejected
			self.selection = 0
		else:
			self.selection = 1
		
	def plot(self,index):
				
		fig = plt.figure()
		ax = fig.add_subplot()
		plt.xlabel("Response Surface estimation")
		plt.ylabel("Direct Simulation")
		plt.plot(np.asarray(self.y_Test_simulation),np.asarray(self.y_Test_Predict),"k.",ms=8,label=f"Testing (max error = {self.ytestMaxError :.3f}%, mean error = {self.yTestMeanError:.3f}%)")
		#plt.scatter(np.asarray(yData_testing), np.asarray(Sim_value_testing), color="none", edgecolor="black")
		plt.scatter(np.asarray(self.Y), np.asarray(self.resFramWrk), color="none", edgecolor="green",label=f"Training (max error = {self.MaxError:.3f}%, mean error = {self.MeanError:.3f}%)")
	
		x = np.linspace(0,500,1000)
		plt.xlim(min(np.asarray(self.Y))*0.98,max(np.asarray(self.Y))*1.02)
		plt.ylim(min(np.asarray(self.resFramWrk))*0.98,max(np.asarray(self.resFramWrk))*1.02)
		plt.plot(x,x,"-",label="parity line")
		plt.legend(loc="upper left")
		print(os.getcwd(),index)
		plt.savefig('../Plots/Parity_plot_case_'+str(self.case_index)+'_training.png',bbox_inches="tight")	
	
	def create_gauss_response_surface(self):
		kernel = DotProduct() + WhiteKernel()
		self.gpr = GaussianProcessRegressor(kernel=kernel,
											random_state=0).fit(self.X, self.Y)
		self.gpr_score = self.gpr.score(self.X,self.Y)
	#def create_neural_network(self):
		self.gpr_yfit = self.gpr.predict(self.X)
	
	def create_Isotonic_response_surface(self):
		self.huber = QuantileRegressor(quantile=0.8).fit(self.X, self.Y)
		self.huber_yfit = self.huber.predict(self.X)
	
	def create_HuberRegressor_response_surface(self):
		self.huber = HuberRegressor().fit(self.X, self.Y)
		self.huber_yfit = self.huber.predict(self.X)
		"""
		self.resFramWrk = []
		for i in self.X:
		    self.resFramWrk.append(self.huber.predict(i))
		
		fileResponse = open(os.getcwd()+'/Data/ResponseSurface/FlaMan_Response_comparison_.csv','w')
		simVSresp  = "FlameMaster,Response Surface,Error_(ln(tau)),Error(Tau)\n"
		self.RMS_error = []
		self.error = []
		self.relative_error = []

		TraError = []
		for i in range(len(self.resFramWrk)):
			self.error.append(abs(self.Y[i]-self.resFramWrk[i]))
			self.RMS_error.append(abs(self.Y[i]-self.resFramWrk[i])**2)
			self.relative_error.append(abs(self.Y[i]-self.resFramWrk[i])/(self.Y[i]))
			TraError.append(1-np.exp(-(self.Y[i]-self.resFramWrk[i])))
			simVSresp +='{},{},{},{}\n'.format(self.Y[i],self.resFramWrk[i],(self.Y[i]-self.resFramWrk[i])/self.Y[i],1-np.exp(-(self.Y[i]-self.resFramWrk[i])))	
		fileResponse.write(simVSresp)		
		fileResponse.close() 
		"""
		"""
		error details
		"""
		"""
		self.MaxError = max(self.error)
		self.MeanError = statistics.mean(self.error)
		
		self.RMS = math.sqrt(sum(self.RMS_error)/len(self.RMS_error))
		self.MaxError = max(self.relative_error)
		self.MeanError = statistics.mean(self.relative_error)
		"""
	
	def get_response_surface(self,xData):
		self.X = xData
		self.BTrsMatrix = self.MatPolyFitTransform()
		self.Q, self.R = np.linalg.qr(self.BTrsMatrix)
		y = np.dot(np.transpose(self.Q),self.Y)
		self.coeff = np.linalg.solve(self.R,np.transpose(y))
		"""
		The simulation data is
		
		self.Y
		"""
		"""
		Writing the response surface for checking the co-efficients
		"""
		rr = open(os.getcwd()+'/Data/ResponseSurface/responsecoef.csv','w')		
		res = "Coefficients\n"
		for i in self.coeff:
			res +='{}\n'.format(float(i))
		rr.write(res)
		rr.close()	
		
		self.resFramWrk = []
		for i in self.X:
		    self.resFramWrk.append(self.evaluate(i))
		
		fileResponse = open(os.getcwd()+'/Data/ResponseSurface/FlaMan_Response_comparison_.csv','w')
		simVSresp  = "FlameMaster,Response Surface,Error_(ln(tau)),Error(Tau)\n"
		self.RMS_error = []
		self.error = []
		self.relative_error = []

		TraError = []
		for i in range(len(self.resFramWrk)):
			self.error.append(abs(self.Y[i]-self.resFramWrk[i]))
			self.RMS_error.append(abs(self.Y[i]-self.resFramWrk[i])**2)
			self.relative_error.append(abs(self.Y[i]-self.resFramWrk[i])/(self.Y[i]))
			TraError.append(1-np.exp(-(self.Y[i]-self.resFramWrk[i])))
			simVSresp +='{},{},{},{}\n'.format(self.Y[i],self.resFramWrk[i],(self.Y[i]-self.resFramWrk[i])/self.Y[i],1-np.exp(-(self.Y[i]-self.resFramWrk[i])))	
		fileResponse.write(simVSresp)		
		fileResponse.close() 
		
		"""
		error details
		"""
		#self.MaxError = max(self.error)
		#self.MeanError = statistics.mean(self.error)
		
		self.RMS = math.sqrt(sum(self.RMS_error)/len(self.RMS_error))
		self.MaxError = max(self.relative_error)
		self.MeanError = statistics.mean(self.relative_error)
		del self.X,self.BTrsMatrix,self.Q, self.R
		
	def create_response_surface(self):
		
		if os.getcwd()+'/Data/ResponseSurface/responsecoef_case-'+str(self.case_index)+'.csv' is True:
			#print("The file is there") 
			f = open(os.getcwd()+'/Data/ResponseSurface/responsecoef_case-'+str(self.case_index)+'.csv','r').readlines()
			self.coeff = np.asarray([float(i) for i in f[1:]])
			#print(self.coeff)
			#print("The file is there")
		
		else:
			self.BTrsMatrix = self.MatPolyFitTransform()
			#print(np.shape(self.X))
			self.Q, self.R = np.linalg.qr(self.BTrsMatrix)
			#print(np.shape(self.Q))
			#print(np.shape(self.R))
			y = np.dot(np.transpose(self.Q),self.Y)
			self.coeff = np.linalg.solve(self.R,np.transpose(y))
			#print(self.coeff)
			"""
			Writing the response surface for checking the co-efficients
			"""
			rr = open(os.getcwd()+'/Data/ResponseSurface/responsecoef_case-'+str(self.case_index)+'.csv','w')		
			res = "Coefficients\n"
			for i in self.coeff:
				res +='{}\n'.format(float(i))
			rr.write(res)
			rr.close()		
			
			if self.order == 2:
				self.zero,self.a,self.b = self.resCoeffTransform(self.order)

			del self.BTrsMatrix,self.Q, self.R
		#raise AssertionError("Fast PRS")
		#coeff_size = np.asarray(self.BTrsMatrix[0]).flatten()
		#guess = np.ones(len(coeff_size))
		#x = minimize(self.objective,guess,method="SLSQP")
		#self.coeff = x.x
		
		
			
			

	
		"""
		factor of evaluating the performance of the
		response surface 
		"""
		"""
		norm_Ax_b_0 = np.linalg.norm(abs(np.dot(self.R,self.coeff)-y),ord = 0)
		norm_b_0 = np.linalg.norm(abs(y),ord = 0)
		norm_Ax_b_1 = np.linalg.norm(abs(np.dot(self.R,self.coeff)-y),ord = 1)
		norm_b_1 = np.linalg.norm(abs(y),ord = 1)
		norm_Ax_b_2 = np.linalg.norm(abs(np.dot(self.R,self.coeff)-y),ord = 2)
		norm_b_2 = np.linalg.norm(abs(y),ord = 2)
		"""
		"""
		Writing the prediction of responseSurface
		"""
		
		self.resFramWrk = []
		for i in self.X:
		    self.resFramWrk.append(self.evaluate(i))
		
		#print(self.resFramWrk)
		fileResponse = open(os.getcwd()+'/Data/ResponseSurface/FlaMan_Response_comparison_.csv','w')
		simVSresp  = "FlameMaster,Response Surface,Error_(ln(tau)),Error(Tau)\n"
		self.RMS_error = []
		self.error = []
		self.relative_error = []

		TraError = []
		for i in range(len(self.resFramWrk)):
			self.error.append(abs(self.Y[i]-self.resFramWrk[i]))
			self.RMS_error.append(abs(self.Y[i]-self.resFramWrk[i])**2)
			self.relative_error.append(abs(self.Y[i]-self.resFramWrk[i])/(self.Y[i]))
			TraError.append(1-np.exp(-(self.Y[i]-self.resFramWrk[i])))
			simVSresp +='{},{},{},{}\n'.format(self.Y[i],self.resFramWrk[i],(self.Y[i]-self.resFramWrk[i])/self.Y[i],1-np.exp(-(self.Y[i]-self.resFramWrk[i])))	
		fileResponse.write(simVSresp)		
		fileResponse.close() 
		
		"""
		error details
		"""
		self.MaxError = max(self.error)
		self.MeanError = statistics.mean(self.error)
		
		self.RMS = math.sqrt(sum(self.RMS_error)/len(self.RMS_error))
		self.MaxError = max(self.relative_error)
		self.MeanError = statistics.mean(self.relative_error)
		del self.X
	
	def MatPolyFitTransform(self):
		BTrsMatrix = []
		for i in self.X:
			tow = self.order
			row = i
			row_ = []
			row_.append(1)
			if tow > 0:		
				for i in row:				
					row_.append(i)
				tow = tow - 1

			if tow > 0:
				for i,j in enumerate(row): 
					for k in row[i:]:					
						row_.append(j*k)
				tow = tow - 1

			if tow > 0:		
				for i in row:
					for j in row[row.index(i):]: 
						for k in row[row.index(j):]:					
							row_.append(i*j*k)
				tow = tow - 1					

			if tow > 0:
				for i,j in enumerate(row):
					for j in row[row.index(i):]: 
						for k in row[row.index(j):]:					
							for l in row[row.index(k):]:
								row_.append(i*j*k*l)
				tow = tow - 1

			if tow > 0:
				for i in row:
					for j in row[row.index(i):]: 
						for k in row[row.index(j):]:					
							for l in row[row.index(k):]:
								for m in row[row.index(l):]:
									row_.append(i*j*k*l*m)
				tow = tow - 1
			BTrsMatrix.append(row_)
		return BTrsMatrix
	
	def evaluate(self,x):
		BZeta = x
		coeff = self.coeff
		tow = self.order
		row = BZeta
		val = coeff[0]
		count = 1
		if tow > 0:		
			for i in BZeta:				
				val+=coeff[count]*i
				count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(BZeta): 
				for k in BZeta[i:]:
					if count < len(coeff):					
						val+=coeff[count]*j*k
						count +=1
			tow = tow - 1

		if tow > 0:		
			for i,j in enumerate(BZeta):
				for k in BZeta[BZeta.index(j):]: 
					for m in row[row.index(k):]:
						if count < len(coeff):					
							val+=coeff[count]*j*k*m
							count +=1
			tow = tow - 1					

		if tow > 0:
			for i,j in enumerate(BZeta):
				for k,l in enumerate(BZeta[i:]): 
					for m,n in enumerate(BZeta[k:]):					
						for o in BZeta[m:]:
							if count < len(coeff):
								val+=coeff[count]*j*l*n*o
								count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(BZeta):
				for k,l in enumerate(BZeta[i:]): 
					for m,n in enumerate(BZeta[k:]):					
						for o,p in enumerate(BZeta[m:]):
							for q in BZeta[o:]:
								if count<len(coeff):
									val+=coeff[count]*q*p*n*l*j
									count +=1
			tow = tow - 1
			
		return val
	
	def evaluateResponse(self,x,cov_x=None):
		val = self.zero
		val += np.dot(self.a,x)
		
		if cov_x is not None:
			a_times_cov = np.dot(self.a,cov_x)
			variance = np.dot(self.a,a_times_cov.T)
        
		if self.b is not None:
			b_times_x = np.asarray(np.dot(self.b,x)).flatten()
			val += np.dot(b_times_x.T,x)
			if cov_x is not None:
				b_times_cov = np.dot(self.b,cov_x)
				variance += 2*np.trace(np.dot(b_times_cov,b_times_cov))
         
		if cov_x is not None:
			computed_unc = math.sqrt(variance)
			return val,computed_unc
		return val
	
	def Jacobian(self,x):
		j = []
		x = list(x)
		#print("\tx{} (type-{}) in response surface is\n".format(x,type(x)))
		for i,opt in enumerate(x):
			j.append(self.jacobian_element(self.coeff,x,opt,i,2))
		return j 
	
	def estimate(self,x):
		val = self.zero
		val += np.dot(self.a,x)
		
		response_grad = self.a
		
		if not(self.b is None):
			b_times_x = np.asarray(np.dot(self.b,x)).flatten()
			val += np.dot(b_times_x.T,x)
			response_grad += 2*b_times_x
            
		return val,response_grad
		
	def resCoeffTransform(self,order):
		#x = self.xdata
		coeff = self.coeff
		tow = order
		row = self.X[0]
		zero = []
		a = []
		b = []
		zero.append(coeff[0])
		count = 1
		if tow > 0:		
			for _ in row:				
				a.append(coeff[count])
				count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(row): 
				temp = []
				for k,l in enumerate(row[i:]):
					if count < len(coeff):				
						temp.append(coeff[count])
						count +=1
				b.append(list(np.zeros(len(row)-len(temp)))+temp)					
			tow = tow - 1
		return float(self.coeff[0]),np.asarray(a),np.matrix(b)
	
	def objective(self,z):
		"""
		(Ax-y) + penalty = 0
		"""
		A = self.BTrsMatrix
		prediction = A.dot(z)
		simulation = self.actualValue
		residual = (prediction - simulation)
		obj = np.dot(residual,residual.T)
		for i in z:
			obj+=i**2
		return obj
		
	def model_response_uncertainty(self,x,cov,order):
			
		coeff = self.resCoef
		tow = order
		val = 0
		a = 0
		b_ii = 0
		b_ij = 0
		row = BZeta
		count = 1
		if tow > 0:		
			for i in BZeta:				
				a += coeff[count]**2
				count +=1
			tow = tow - 1

		if tow > 0:
			for i,j in enumerate(BZeta): 
				for k,l in enumerate(BZeta[i:]):
					if count < len(coeff):
						if i == k:
							b_ii+=coeff[count]**2
						else:
							b_ij+=coeff[count]**2					
						count +=1
			tow = tow - 1
		
		#print(a)
		#print(b_ii)
		#print(b_ij)
		val = a + 2*b_ii+b_ij
		#val = np.linalg.norm(np.dot(np.asarray(a),cov))**2+ np.linalg.norm(np.dot(cov.T,np.dot(np.matrix(b),cov)),'fro')**2
		### Calculating the variance for each response surface
		#print("Len of x mdres is {}".format(len(self.resCoef)))
		#print("Count is {}\n".format(count))	
		return np.sqrt(val)
	
	
	def jacobian_element(self,coeff,BZeta,x,ind,resp_order):
		"""
		J = a + bx
		"""
		tow = resp_order
		row = BZeta   
		val = 0
		count = 1
		index=[]
		"""
		find a
		"""
		if tow > 0:
			for i,j in enumerate(BZeta):
				
				if i == ind:
					val+=coeff[count]
					index.append(count)
					count +=1
				
				else:
					count +=1
			tow = tow - 1

		"""
		find b*x
		"""
		if tow > 0:
			for i,j in enumerate(BZeta):
				l = i 
				for k in BZeta[i:]:
					#if count < len(coeff):
					if i == ind and l == ind:
						val += 2*coeff[count]*k
						index.append(count)
						#print("{}{}\t{}\n".format(i,l,count))
						count+=1
						l+=1
					elif i == ind and l!=ind:
						val += coeff[count]*k
						index.append(count)
						#print("{}{}\t{}\n".format(i,BZeta.index(k),count))
						count+=1
						l+=1
					elif i!= ind and l==ind:
						val += coeff[count]*j
		#				print("{}{}\t{}\n".format(i,BZeta.index(k),count))
						index.append(count)
						count+=1
						l+=1
					else:
						#print("uncounted {}{}\t{}\n".format(i,BZeta.index(k),count))
						count+=1
						l+=1
			tow = tow - 1
		#print(val)
		#print("\t\tThe index for {} in x[{}]\n{}\n".format(x,ind,index))
		return val
	
