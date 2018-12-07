
.. image:: https://img.shields.io/github/forks/badges/shields.svg?style=social&label=Fork   
    :target: https://github.com/DataXujing/LN0SCIs/

.. image:: https://img.shields.io/pypi/pyversions/Django.svg   
	:target: https://pypi.python.org/pypi/LN0SCIs

.. image:: https://img.shields.io/cran/v/devtools.svg   
	:target: https://CRAN.R-project.org/package=LN0SCIs

.. image:: https://img.shields.io/pypi/v/nine.svg   
	:target: https://pypi.python.org/pypi/LN0SCIs

.. image:: https://github.com/DataXujing/LN0SCIs/raw/master/pic/log.png
    :align: right

LN0SCIs
===============

**Jing Xu, Xinmin Li, Hua Liang**



Introduction
---------------

This Python package based on the paper of Simultaneous Confidence Intervals for Ratios of Means of Log-normal Populations with Zeros by Xu et al. It provides some methods for construct simultaneous confidence intervals for ratios of means of Log-normal populations with excess zeros. At last, we select 4 excellent methods which based on generalized pivotal quantity with order statistics and two-step MOVER intervals. For the convenience of use, we make a Python package called LN0SCIs, and it also has a R version package on CRAN: https://CRAN.R-project.org/package=LN0SCIs

+ If you are a R User, you can install in your R kernal by Github:

  - `devtools::install_github('DataXujing/LN0SCIs')`

+ Or you can also install by CRAN:

  - `install.packages('LN0SCIs')`

+ If you are a Python user, you can 

  - `pip install LN0SCIs`



Methods
------------

We provaide four main functions in our LN0SCIs packages, FGW(),FGH(),MOVERW() and MOVERH(), if you want to deep understanding these four methods, you can read our paper: *Simultaneous Confidence Intervals for Ratios of Means of Log-normal Populations with Zeros*. the code we trust in GitHub. If you want to know how to realize them, you can read the source code.


Examples
---------

+ FGW()

::


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

  
::

	====================Method: FGW=====================
	The Simultaneous Confidence Intervals are:          
             The1th CIs            The2th CIs            The3th CIs
	0  【-0.843638,0.789044】  【-0.629208,1.075959】  【-0.604469,1.158544】
	**********************Time**************************
	The cost time is:0 secs
	====================Method: FGW=====================
	The Simultaneous Confidence Intervals are:          
             The1th CIs           The2th CIs           The3th CIs  \
	0  【-0.912169,1.578679】  【-1.02404,0.812882】  【-0.83778,1.382352】   

             The4th CIs            The5th CIs           The6th CIs  
	0  【-1.597962,0.650222】  【-1.337939,1.203199】  【-0.546039,1.25945】  
	**********************Time**************************
	The cost time is:0 secs


+ FGH()

::

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

::

	====================Method: FGH=====================
	The Simultaneous Confidence Intervals are:          
             The1th CIs            The2th CIs            The3th CIs
	0  【-0.992276,1.455247】  【-0.703231,1.372774】  【-1.005873,1.124758】
	**********************Time**************************
	The cost time is:0 secs
	====================Method: FGH=====================
	The Simultaneous Confidence Intervals are:          
            The1th CIs            The2th CIs            The3th CIs  \
	0  【-1.62426,0.624984】  【-1.514528,0.553936】  【-1.565943,0.911157】   

            The4th CIs            The5th CIs           The6th CIs  
	0  【-0.66646,1.010746】  【-0.829753,1.269381】  【-0.762683,1.07889】  
	**********************Time**************************
	The cost time is:0 secs


+ MOVERW()


::


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


::


	====================Method: FGH=====================
	The Simultaneous Confidence Intervals are:          
             The1th CIs            The2th CIs            The3th CIs
	0  【-1.103496,1.211033】  【-1.030952,0.888781】  【-1.314926,1.059975】
	**********************Time**************************
	The cost time is:0 secs
	====================Method: FGH=====================
	The Simultaneous Confidence Intervals are:          
            The1th CIs            The2th CIs            The3th CIs  \
	0  【-1.68825,0.349316】  【-1.270833,1.236153】  【-1.304731,1.053776】   

             The4th CIs            The5th CIs            The6th CIs  
	0  【-0.349427,1.679719】  【-0.364992,1.484843】  【-1.294225,1.071433】  
	**********************Time**************************
	The cost time is:0 secs


+ MOVERH()


::


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


::

	====================Method: FGH=====================
	The Simultaneous Confidence Intervals are:          
             The1th CIs            The2th CIs          The3th CIs
	0  【-1.013305,0.765726】  【-1.152934,0.823283】  【-0.914194,0.8239】
	**********************Time**************************
	The cost time is:0 secs
	====================Method: FGH=====================
	The Simultaneous Confidence Intervals are:          
             The1th CIs            The2th CIs           The3th CIs  \
	0  【-0.681666,1.693927】  【-0.750657,1.458978】  【-1.21012,0.855608】   

             The4th CIs            The5th CIs            The6th CIs  
	0  【-1.302431,1.003355】  【-1.762379,0.407925】  【-1.527028,0.467458】  
	**********************Time**************************
	The cost time is:0 secs	





Supports
-----------

Tested on Python 2.7, 3.5, 3.6

* pip install LN0SCIs
* Download: https://pypi.python.org/pypi/LN0SCIs
* Documentation: https://github.com/DataXujing/LN0SCIs
* It has a R packages version which we have created, details you can see:  https://CRAN.R-project.org/package=LN0SCIs

you can log in Xujing's home page: https://dataxujing.coding.me or https://dataxujing.github.io to find the author(s), and if you want to learn more about simultaneous confidence intervals for the mixture distribution, you shou read the paper: Simulataneous Confidence Intervals for ratios of Means of Log-normal Populations with Zeros, which written by Jing Xu, Xinmin Li, and Hua Liang.

