import codecs
import os
import sys

try:
	from setuptools import setup
except:
	from distutils.core import setup



def read(fname):
	return codecs.open(os.path.join(os.path.dirname(__file__), fname),encoding='utf-8').read()

long_des = read("README.rst")
    
platforms = ['Linux/Windows']
classifiers = [
    'Development Status :: 3 - Alpha',
    'Topic :: Text Processing',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
]

install_requires = [
    'numpy>=1.13.0',
    'pandas>=0.21.0',
    'datetime'
]

    
setup(name='LN0SCIs',
      version='0.1.2',
      description='Simultaneous CIs for Ratios of Means of Log-Normal Populations with Zeros',
      long_description=long_des,
      py_modules=['LN0SCIs'],
      author = "Jing Xu, Xinmin Li, Hua Liang",  
      author_email = "274762204@qq.com" ,
      url = "https://dataxujing.github.io/" ,
      license="Apache License, Version 2.0",
      platforms=platforms,
      classifiers=classifiers,
      install_requires=install_requires
      
      )   
      
      
