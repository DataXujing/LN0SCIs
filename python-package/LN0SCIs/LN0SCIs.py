# -*- coding: utf-8 -*-

import numpy as np
import datetime
import pandas as pd


#basic Fuc
def rndMixture(N0,N1,p,N):
	u = np.random.uniform(low=0.0, high=1.0, size=N)
	r = list()
	for i in range(N):
		if u[i] < p:
			r.append(np.random.beta(N0+1,N1))
		else:
			r.append(np.random.beta(N0,N1+1))

	return r



#FGW-methods
def FGW(n,p,mu,sigma,N,C2=np.array([[-1,1,0],[-1,0,1],[0,-1,1]]) ,alpha=0.05):
	'''
	====================================================================
	The SCIs Based on FGW Methods
	Author(s): Jing Xu, Xinmin Li, Hua Liang
	=====================================================================
	A method based on generalized pivotal quantity with order statistics
	(also see help(FGH)) to construct the simultaneous confidence intervals 
	for Ratios of Means of Log-normal Populations with excess Zeros.
	====================================================================

	Useage:

		FGW(n,p,mu,sigma,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1)),alpha=0.05)

	Params:

	n     The sample size of the mixture distributions,must be an integer vector.
	p     The zero probability of the mixture distribution,it has the same length 
	      to the n params.
	mu    The mean of the non-zero samples,which after log-transformation.
	sigma The variance of the non-zero samples,which after log-transformation.
	N     The number of independent generated data sets.
	C2    Matrix C,You can refer to the paper of Xu et al. for specific forms.
	alpha The confidence level,it always set alpha=0.5

	Return:

	    The method will return the Simultaneous Confidence Intervals(SCIs) and the 
		 time consuming.

	Examples:

	    import numpy as np
		from LN0SCIs import *
		#Example1:
		alpha = 0.05
		p = np.array([0.2,0.2,0.2])
		n = np.array([30,30,30])
		mu = np.array([0,0,0])
		sigma = np.array([1,1,1])
		N = 1000
		FGW(n,p,mu,sigma,N)
		#Example2:
		p = np.array([0.1,0.1,0.1,0.1])
		n = np.array([30,30,30,30])
		mu = np.array([0,0,0,0])
		sigma = np.array([1,1,1,1])
		C2 = np.array([[-1,1,0,0],[-1,0,1,0],[-1,0,0,1],[0,-1,1,0],[0,-1,0,1],[0,0,-1,1]])
		N = 1000
		FGW(n,p,mu,sigma,N,C2 = C2)

	Note: 
		We also create a R package do the them thing: https://CRAN.R-project.org/package=LN0SCIs


	'''
	t1 = datetime.datetime.now()

	nk0 = dict()
	nk1 = dict()
	ykbar = dict()
	sksq = dict()
	Zk = dict()
	Uk = dict()
	ZkW = dict()
	for i in range(len(n)):
		nk0['n'+str(i)+'0'] = np.random.binomial(n[i],p[i])
		nk1['n'+str(i)+'1'] = n[i] - nk0['n'+str(i)+'0']
		ykbar['y'+str(i)+'bar'] = np.random.normal(loc=mu[i],scale=np.float(sigma[i])/np.sqrt(nk1['n'+str(i)+'1']))
		sksq['s'+str(i)+'sq']  = (sigma[i])**2*np.random.chisquare(df=(nk1['n'+str(i)+'1'] - 1))

		Zk['Z'+str(i)] = np.random.normal(loc=0,scale=1,size=N)
		Uk['U'+str(i)] = np.random.chisquare(df=(nk1['n'+str(i)+'1'] - 1),size=N)
		ZkW['Z'+str(i)+'W'] = np.random.normal(loc=0,scale=1,size=N)


	TPW = dict()
	TW = dict()
	for i in range(len(n)):
		TPW['TPW'+str(i)]  = list()
		TW['TW'+str(i)] = list()

	for k in range(N):
		for i in range(len(n)):
			TPW['TPW'+str(i)].append((nk0['n'+str(i)+'0']+ 0.5*(ZkW['Z'+str(i)+'W'][k])**2)/(n[i]+(ZkW['Z'+str(i)+'W'][k])**2)-\
				(ZkW['Z'+str(i)+'W'][k]*np.sqrt(nk0['n'+str(i)+'0']*(1-np.float(nk0['n'+str(i)+'0'])/n[i])+(ZkW['Z'+str(i)+'W'][k])**2/4))/\
				(n[i]+(ZkW['Z'+str(i)+'W'][k])**2))
			TW['TW'+str(i)].append(np.log(1-TPW['TPW'+str(i)][k])+\
				(ykbar['y'+str(i)+'bar']-Zk['Z'+str(i)][k]*np.sqrt(sksq['s'+str(i)+'sq'])/(np.sqrt(nk1['n'+str(i)+'1']*Uk['U'+str(i)][k]))+\
					sksq['s'+str(i)+'sq']/(2*Uk['U'+str(i)][k])))

 

	TWdf = pd.DataFrame(TW)
	TWarray = np.array(TWdf)
	TWmat = np.mat(TWarray)  #N*3
	r = np.mat(np.zeros((C2.shape[0],N)))
	r1 = np.mat(np.zeros((C2.shape[0],N)))
	r2 = np.mat(np.zeros((C2.shape[0],N)))
	for k in range(N):
		r[:,k] = np.mat(C2)*TWmat[k].T 


	for j in range(C2.shape[0]):
		r1[j,:] = np.sort(r[j,:])
		r2[j,:] = np.array(pd.Series(np.array(r[j,:])[0]).rank())



	mik1 = list()
	mak1 = list()

	for k in range(N):
		mik1.append(np.min(r2[:,k]))
		mak1.append(np.max(r2[:,k]))


	mik2 = np.sort(mik1)
	mak2 = np.sort(mak1)

	kl = np.int(mik2[np.int(N*alpha/2)])
	ku = np.int(mak2[np.int(N*(1-alpha/2))])

	CL = dict()
	for j in range(C2.shape[0]):
		CL['The'+str(j+1)+'th'+' CIs'] = ['【'+str(round(r1[j,kl],6))+','+str(round(r1[j,ku],6))+'】']


	CLdf = pd.DataFrame(CL)
	t2 = datetime.datetime.now()


	print("====================Method: FGW=====================")
	print("The Simultaneous Confidence Intervals are:          ")
	print(CLdf)
	print('**********************Time**************************')
	print('The cost time is:'+str((t2-t1).seconds)+' secs')

#FGH

def FGH(n,p,mu,sigma,N,C2=np.array([[-1,1,0],[-1,0,1],[0,-1,1]]) ,alpha=0.05):
	'''
	====================================================================
	The SCIs Based on FGH Methods
	Author(s): Jing Xu, Xinmin Li, Hua Liang
	=====================================================================
	A method based on generalized pivotal quantity with order statistics
	(also see help(FGW)) to construct the simultaneous confidence intervals 
	for Ratios of Means of Log-normal Populations with excess Zeros.
	====================================================================

	Useage:

		FGH(n,p,mu,sigma,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1)),alpha=0.05)

	Params:

	n     The sample size of the mixture distributions,must be an integer vector.
	p     The zero probability of the mixture distribution,it has the same length 
	      to the n params.
	mu    The mean of the non-zero samples,which after log-transformation.
	sigma The variance of the non-zero samples,which after log-transformation.
	N     The number of independent generated data sets.
	C2    Matrix C,You can refer to the paper of Xu et al. for specific forms.
	alpha The confidence level,it always set alpha=0.5

	Return:

	    The method will return the Simultaneous Confidence Intervals(SCIs) and the 
		 time consuming.

	Examples:

	    import numpy as np
		from LN0SCIs import *
		#Example1:
		alpha = 0.05
		p = np.array([0.2,0.2,0.2])
		n = np.array([30,30,30])
		mu = np.array([0,0,0])
		sigma = np.array([1,1,1])
		N = 1000
		FGH(n,p,mu,sigma,N)
		#Example2:
		p = np.array([0.1,0.1,0.1,0.1])
		n = np.array([30,30,30,30])
		mu = np.array([0,0,0,0])
		sigma = np.array([1,1,1,1])
		C2 = np.array([[-1,1,0,0],[-1,0,1,0],[-1,0,0,1],[0,-1,1,0],[0,-1,0,1],[0,0,-1,1]])
		N = 1000
		FGH(n,p,mu,sigma,N,C2 = C2)

	Note: 
		We also create a R package do the them thing: https://CRAN.R-project.org/package=LN0SCIs


	'''

	t1 = datetime.datetime.now()
	nk0 = dict()
	nk1 = dict()
	ykbar = dict()
	sksq = dict()
	Zk = dict()
	Uk = dict()
	TPH = dict()
	for i in range(len(n)):
		nk0['n'+str(i)+'0'] = np.random.binomial(n[i],p[i])
		nk1['n'+str(i)+'1'] = n[i] - nk0['n'+str(i)+'0']
		ykbar['y'+str(i)+'bar'] = np.random.normal(loc=mu[i],scale=np.float(sigma[i])/np.sqrt(nk1['n'+str(i)+'1']))
		sksq['s'+str(i)+'sq']  = (sigma[i])**2*np.random.chisquare(df=(nk1['n'+str(i)+'1'] - 1))

		Zk['Z'+str(i)] = np.random.normal(loc=0,scale=1,size=N)
		Uk['U'+str(i)] = np.random.chisquare(df=(nk1['n'+str(i)+'1'] - 1),size=N)
		TPH['TPH'+str(i)] = rndMixture(nk0['n'+str(i)+'0'],nk1['n'+str(i)+'1'],0.5,N)

	TW = dict()
	for i in range(len(n)):
		TW['TW'+str(i)] = list()

	for k in range(N):
		for i in range(len(n)):
			TW['TW'+str(i)].append(np.log(1-TPH['TPH'+str(i)][k])+\
				(ykbar['y'+str(i)+'bar']-Zk['Z'+str(i)][k]*np.sqrt(sksq['s'+str(i)+'sq'])/\
					(np.sqrt(nk1['n'+str(i)+'1']*Uk["U"+str(i)][k])) +\
					sksq['s'+str(i)+'sq']/(2*Uk['U'+str(i)][k])))


	TWdf = pd.DataFrame(TW)
	TWarray = np.array(TWdf)
	TWmat = np.mat(TWarray)  #N*3
	r = np.mat(np.zeros((C2.shape[0],N)))
	r1 = np.mat(np.zeros((C2.shape[0],N)))
	r2 = np.mat(np.zeros((C2.shape[0],N)))
	for k in range(N):
		r[:,k] = np.mat(C2)*TWmat[k].T 


	for j in range(C2.shape[0]):
		r1[j,:] = np.sort(r[j,:])
		r2[j,:] = np.array(pd.Series(np.array(r[j,:])[0]).rank())



	mik1 = list()
	mak1 = list()

	for k in range(N):
		mik1.append(np.min(r2[:,k]))
		mak1.append(np.max(r2[:,k]))


	mik2 = np.sort(mik1)
	mak2 = np.sort(mak1)

	kl = np.int(mik2[np.int(N*alpha/2)])
	ku = np.int(mak2[np.int(N*(1-alpha/2))])

	CL = dict()
	for j in range(C2.shape[0]):
		CL['The'+str(j+1)+'th'+' CIs'] = ['【'+str(round(r1[j,kl],6))+','+str(round(r1[j,ku],6))+'】']


	CLdf = pd.DataFrame(CL)
	t2 = datetime.datetime.now()


	print("====================Method: FGH=====================")
	print("The Simultaneous Confidence Intervals are:          ")
	print(CLdf)
	print('**********************Time**************************')
	print('The cost time is:'+str((t2-t1).seconds)+' secs')


#MOVERW

def MOVERW(n,p,mu,sigma,N,C2=np.array([[-1,1,0],[-1,0,1],[0,-1,1]]) ,alpha=0.05):
	'''
	====================================================================
	The SCIs Based on MOVERW Methods
	Author(s): Jing Xu, Xinmin Li, Hua Liang
	=====================================================================
	A method based on two-step MOVER intervals(also see MOVERH)
	to construct the simultaneous confidence intervals for Ratios of Means of 
	Log-normal Populations with Zeros.
	====================================================================

	Useage:

		MOVERW(n,p,mu,sigma,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1)),alpha=0.05)

	Params:

	n     The sample size of the mixture distributions,must be an integer vector.
	p     The zero probability of the mixture distribution,it has the same length 
	      to the n params.
	mu    The mean of the non-zero samples,which after log-transformation.
	sigma The variance of the non-zero samples,which after log-transformation.
	N     The number of independent generated data sets.
	C2    Matrix C,You can refer to the paper of Xu et al. for specific forms.
	alpha The confidence level,it always set alpha=0.5

	Return:

	    The method will return the Simultaneous Confidence Intervals(SCIs) and the 
		time consuming.

	Examples:

	    import numpy as np
		from LN0SCIs import *
		#Example1:
		alpha = 0.05
		p = np.array([0.2,0.2,0.2])
		n = np.array([30,30,30])
		mu = np.array([0,0,0])
		sigma = np.array([1,1,1])
		N = 1000
		MOVERW(n,p,mu,sigma,N)
		#Example2:
		p = np.array([0.1,0.1,0.1,0.1])
		n = np.array([30,30,30,30])
		mu = np.array([0,0,0,0])
		sigma = np.array([1,1,1,1])
		C2 = np.array([[-1,1,0,0],[-1,0,1,0],[-1,0,0,1],[0,-1,1,0],[0,-1,0,1],[0,0,-1,1]])
		N = 1000
		MOVERW(n,p,mu,sigma,N,C2 = C2)

	Note: 
		We also create a R package do the them thing: https://CRAN.R-project.org/package=LN0SCIs


	'''

	t1 = datetime.datetime.now()

	nk0 = dict()
	nk1 = dict()
	pkh = dict()
	ykbar = dict()
	sksq = dict()
	sksq1 = dict()
	Zk = dict()
	Uk = dict()
	ZkW = dict()
	for i in range(len(n)):
		nk0['n'+str(i)+'0'] = np.random.binomial(n[i],p[i])
		nk1['n'+str(i)+'1'] = n[i] - nk0['n'+str(i)+'0']
		pkh['p'+str(i)+'h'] = np.float(nk0['n'+str(i)+'0'])/n[i]

		ykbar['y'+str(i)+'bar'] = np.random.normal(loc=mu[i],scale=np.float(sigma[i])/np.sqrt(nk1['n'+str(i)+'1']))
		sksq['s'+str(i)+'sq']  = (sigma[i])**2*np.random.chisquare(df=(nk1['n'+str(i)+'1'] - 1))
		sksq1['s'+str(i)+'sq1'] = np.float(sksq['s'+str(i)+'sq'])/(nk1['n'+str(i)+'1']-1)

		Zk['Z'+str(i)] = np.random.normal(loc=0.0,scale=1.0,size=N)
		Uk['U'+str(i)] = np.random.chisquare(df=(nk1['n'+str(i)+'1'] - 1),size=N)
		ZkW['Z'+str(i)+'W'] = np.random.normal(loc=0.0,scale=1.0,size=N)

	TPW = dict()
	TW = dict()

	for i in range(len(n)):
		TPW['TPW'+str(i)] = list()
		TW['TW'+str(i)] = list()

	for k in range(N):
		for i in range(len(n)):
			TPW['TPW'+str(i)].append((nk0['n'+str(i)+'0']+0.5*(ZkW['Z'+str(i)+'W'][k])**2)/\
				(n[i]+(ZkW['Z'+str(i)+'W'][k])**2)-\
				(ZkW['Z'+str(i)+'W'][k]*np.sqrt(nk0['n'+str(i)+'0']*(1-np.float(nk0['n'+str(i)+'0'])/n[1]) +\
					(ZkW['Z'+str(i)+'W'][k])**2/4))/\
				(n[i]+(ZkW['Z'+str(i)+'W'][k])**2))

			TW['TW'+str(i)].append( np.log(1-TPW['TPW'+str(i)][k])+\
				(ykbar['y'+str(i)+'bar']-Zk['Z'+str(i)][k]*np.sqrt(sksq['s'+str(i)+'sq'])/\
					(np.sqrt(nk1['n'+str(i)+'1']*Uk['U'+str(i)][k]))+\
					sksq['s'+str(i)+'sq']/(2*Uk['U'+str(i)][k])))


        

	TWs = dict()
	for i in range(len(n)):
		TWs['TWs'+str(i)] = np.sort(TW['TW'+str(i)])

	Llimit = list()
	Ulimit = list()
	et = list()
	for i in range(len(n)):
		Llimit.append(TWs['TWs'+str(i)][np.int(N*alpha/(2*C2.shape[0]))])
		Ulimit.append(TWs['TWs'+str(i)][np.int(N*(1-alpha/(2*C2.shape[0])))])
		et.append(np.log(1-pkh['p'+str(i)+'h'])+ykbar['y'+str(i)+'bar']+0.5*sksq1['s'+str(i)+'sq1'])

	Ades = np.mat(C2)*np.mat(et).T
	CLs = dict()
	for j in range(C2.shape[0]):
		CLs['The'+str(j+1)+'th CIs'] = ['【'+ str(round((np.float(Ades[j]) - \
			np.sqrt((et[np.argwhere(C2[j,:]==1)[0][0]] -\
				Llimit[np.argwhere(C2[j,:]==1)[0][0]])**2 + (Ulimit[np.argwhere(C2[j,:]==-1)[0][0]]-\
				et[np.argwhere(C2[j,:]==-1)[0][0]])**2)),6)) +','+str(round((Ades[j] +\
		np.sqrt((Ulimit[np.argwhere(C2[j,:]==1)[0][0]]-et[np.argwhere(C2[j,:]==1)[0][0]])**2 +\
			(et[np.argwhere(C2[j,:]==-1)[0][0]]-Llimit[np.argwhere(C2[j,:]==-1)[0][0]])**2)),6)) + '】']


	CLdf = pd.DataFrame(CLs)
	t2 = datetime.datetime.now()


	print("====================Method: FGH=====================")
	print("The Simultaneous Confidence Intervals are:          ")
	print(CLdf)
	print('**********************Time**************************')
	print('The cost time is:'+str((t2-t1).seconds)+' secs')



#MOVERH


def MOVERH(n,p,mu,sigma,N,C2=np.array([[-1,1,0],[-1,0,1],[0,-1,1]]) ,alpha=0.05):
	'''
	====================================================================
	The SCIs Based on MOVERH Methods
	Author(s): Jing Xu, Xinmin Li, Hua Liang
	=====================================================================
	A method based on two-step MOVER intervals(also see MOVERW)
	to construct the simultaneous confidence intervals for Ratios of Means of 
	Log-normal Populations with Zeros.
	====================================================================

	Useage:

		MOVERH(n,p,mu,sigma,N,C2=rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1)),alpha=0.05)

	Params:

	n     The sample size of the mixture distributions,must be an integer vector.
	p     The zero probability of the mixture distribution,it has the same length 
	      to the n params.
	mu    The mean of the non-zero samples,which after log-transformation.
	sigma The variance of the non-zero samples,which after log-transformation.
	N     The number of independent generated data sets.
	C2    Matrix C,You can refer to the paper of Xu et al. for specific forms.
	alpha The confidence level,it always set alpha=0.5

	Return:

	    The method will return the Simultaneous Confidence Intervals(SCIs) and the 
		time consuming.

	Examples:

	    import numpy as np
		from LN0SCIs import *
		#Example1:
		alpha = 0.05
		p = np.array([0.2,0.2,0.2])
		n = np.array([30,30,30])
		mu = np.array([0,0,0])
		sigma = np.array([1,1,1])
		N = 1000
		MOVERH(n,p,mu,sigma,N)
		#Example2:
		p = np.array([0.1,0.1,0.1,0.1])
		n = np.array([30,30,30,30])
		mu = np.array([0,0,0,0])
		sigma = np.array([1,1,1,1])
		C2 = np.array([[-1,1,0,0],[-1,0,1,0],[-1,0,0,1],[0,-1,1,0],[0,-1,0,1],[0,0,-1,1]])
		N = 1000
		MOVERH(n,p,mu,sigma,N,C2 = C2)

	Note: 
		We also create a R package do the them thing: https://CRAN.R-project.org/package=LN0SCIs


	'''

	t1 = datetime.datetime.now()

	nk0 = dict()
	nk1 = dict()
	pkh = dict()
	ykbar = dict()
	sksq = dict()
	sksq1 = dict()
	Zk = dict()
	Uk = dict()
	TPH = dict()

	for i in range(len(n)):
		nk0['n'+str(i)+'0'] = np.random.binomial(n[i],p[i])
		nk1['n'+str(i)+'1'] = n[i] - nk0['n'+str(i)+'0']
		pkh['p'+str(i)+'h'] = np.float(nk0['n'+str(i)+'0'])/n[i]

		ykbar['y'+str(i)+'bar'] = np.random.normal(loc=mu[i],scale=np.float(sigma[i])/np.sqrt(nk1['n'+str(i)+'1']))
		sksq['s'+str(i)+'sq']  = (sigma[i])**2*np.random.chisquare(df=(nk1['n'+str(i)+'1'] - 1))
		sksq1['s'+str(i)+'sq1'] = np.float(sksq['s'+str(i)+'sq'])/(nk1['n'+str(i)+'1']-1)

		Zk['Z'+str(i)] = np.random.normal(loc=0.0,scale=1.0,size=N)
		Uk['U'+str(i)] = np.random.chisquare(df=(nk1['n'+str(i)+'1'] - 1),size=N)
		TPH['TPH'+str(i)] = rndMixture(nk0['n'+str(i)+'0'],nk1['n'+str(i)+'1'],p[i],N)


	TW = dict()
	for i in range(len(n)):
		TW['TW'+str(i)] = list()

	for k in range(N):
		for i in range(len(n)):
			TW['TW'+str(i)].append(np.log(1-TPH['TPH'+str(i)][k])+\
				(ykbar['y'+str(i)+'bar']-Zk['Z'+str(i)][k]*np.sqrt(sksq['s'+str(i)+'sq'])/\
					(np.sqrt(nk1['n'+str(i)+'1']*Uk['U'+str(i)][k])) + sksq['s'+str(i)+'sq']/(2*Uk['U'+str(i)][k])))


	TWs = dict()
	for i in range(len(n)):
		TWs['TWs'+str(i)] = np.sort(TW['TW'+str(i)])


	Llimit = list()
	Ulimit = list()
	et = list()
	for i in range(len(n)):
		Llimit.append(TWs['TWs'+str(i)][np.int(N*alpha/(2*C2.shape[0]))])
		Ulimit.append(TWs['TWs'+str(i)][np.int(N*(1-alpha/(2*C2.shape[0])))])
		et.append(np.log(1-pkh['p'+str(i)+'h'])+ykbar['y'+str(i)+'bar']+0.5*sksq1['s'+str(i)+'sq1'])

	Ades = np.mat(C2)*np.mat(et).T
	CLs = dict()
	for j in range(C2.shape[0]):
		CLs['The'+str(j+1)+'th CIs'] = ['【'+ str(round((np.float(Ades[j]) - \
			np.sqrt((et[np.argwhere(C2[j,:]==1)[0][0]] -\
				Llimit[np.argwhere(C2[j,:]==1)[0][0]])**2 + (Ulimit[np.argwhere(C2[j,:]==-1)[0][0]]-\
				et[np.argwhere(C2[j,:]==-1)[0][0]])**2)),6)) +','+str(round((Ades[j] +\
		np.sqrt((Ulimit[np.argwhere(C2[j,:]==1)[0][0]]-et[np.argwhere(C2[j,:]==1)[0][0]])**2 +\
			(et[np.argwhere(C2[j,:]==-1)[0][0]]-Llimit[np.argwhere(C2[j,:]==-1)[0][0]])**2)),6)) + '】']


	CLdf = pd.DataFrame(CLs)
	t2 = datetime.datetime.now()


	print("====================Method: FGH=====================")
	print("The Simultaneous Confidence Intervals are:          ")
	print(CLdf)
	print('**********************Time**************************')
	print('The cost time is:'+str((t2-t1).seconds)+' secs')




