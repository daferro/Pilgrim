#!/usr/bin/python2.7
'''
*----------------------------------*
| Package    :  common             |
| Module     :  cag                |
| Last Update:  2018/10/06 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the function
for the calculation of the CAG
correction factor
''' 

#===================================#
from physcons import KB
import numpy as np
#===================================#


#===================================#
def calc_cag(ltemp,VadiSpl,lscvt=None):
    # calculate dE
    sAG, VAG =  VadiSpl.get_max()
    if lscvt is None: Ediff = [VAG-VadiSpl(0.0) for T   in ltemp]
    else            : Ediff = [VAG-VadiSpl(s_i) for s_i in lscvt]
    # get cag
    cag = [np.exp(-dE/KB/T) for T,dE in zip(ltemp,Ediff)]
    return Ediff , cag
#===================================#


