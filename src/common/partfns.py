#!/usr/bin/python2.7
'''
*----------------------------------*
| Package    :  common             |
| Module     :  partfns            |
| Last Update:  2018/10/06 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains some functions related
to the calculation of basic partition functions

Functions (partfns.py):
  * pf_partinbox(mass,T)
  * pf_rigidrotor(imoments,T,rotsymnum=1)
  * pf_harmosc1D(angfreq,T,imag=1E10)
  * pf_harmosc(angfreqs,T,imag=1E10)
  * pf_electr(eslist,T)
'''

import numpy    as     np
from   physcons import TWOPI, KB, H, HBAR, PI, ML

#===============================================#
# Partition functions                           #
#===============================================#
def pf_partinbox(mass,T):
    return (TWOPI*KB*mass*T)**(3./2.)/(H**3)
#-----------------------------------------------#
def pf_rigidrotor(imoments,T,rotsymnum=1):
    beta      = KB * T
    rot_const = [(HBAR**2)/Ii/2 for Ii in imoments]
    # linear case
    if len(imoments) == 1: qrr = beta / rot_const[0]
    else                 : qrr = np.sqrt(PI) * np.sqrt(beta**3 / np.prod(rot_const))
    return qrr/rotsymnum
#-----------------------------------------------#
def pf_harmosc1D(angfreq,T,imag=1E10):
    if angfreq < 0: return imag
    exp  = np.exp(-HBAR*angfreq/KB/T)
    qHO  = 1.0/(1.0-exp)
    return qHO
#-----------------------------------------------#
def pf_harmosc(angfreqs,T,imag=1E10):
    qHO = np.prod([pf_harmosc1D(angfreq,T,imag) for angfreq in angfreqs])
    return qHO
#-----------------------------------------------#
def pf_electr(eslist,T):
    pf, beta = 0.0, 1.0 / KB / T
    for deg, relE in eslist:
        pf += deg * np.exp(-relE * beta)
    return pf
#===============================================#

#==========================================================#
#        Calculation of equilibrium/rate constants         #
#==========================================================#
def Qs2Kc(ltemp,QA,QB,VA,VB):
    '''
    Qi; partition function per unit volume
    Kc = [B]/[A]
    '''
    return [QB[idx]/QA[idx] * np.exp(-(VB-VA)/KB/T) for idx,T in enumerate(ltemp)]
#----------------------------------------------------------#
def Kc2GFE(ltemp,Kc,dn=0):
    c0 = 1.0 / (1.0 / ML) # c0 = 1 particle/mL
    return [ -KB * T * np.log(Kc[idx] * (c0**dn)) for idx,T in enumerate(ltemp)]
#----------------------------------------------------------#
def Kc2rate(ltemp,Kc):
    return [(KB*T/H)*Kc[idx] for idx,T in enumerate(ltemp)]
#----------------------------------------------------------#
def rate2Kc(ltemp,k):
    return [k[idx]*H/(KB*T) for idx,T in enumerate(ltemp)]
#----------------------------------------------------------#
def rate2GFE(ltemp,rates,dn=0):
    Kc  = rate2Kc(ltemp,rates)
    GFE = Kc2GFE(ltemp,Kc,dn)
    return GFE
#==========================================================#

if __name__ == '__main__': 
    import fncs
    fncs.print_funcs_in_module()

