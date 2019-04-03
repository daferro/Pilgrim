#!/usr/bin/python2.7
'''
*----------------------------------*
| Package    :  common             |
| Module     :  steepdesc          |
| Last Update:  2018/10/06 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different functions
related to the calculation of a steepest
descent path
'''

#=======================================================#
import steepdesc  as sd
import fncs       as fncs
import Exceptions as Exc
from   Molecule   import Molecule
from   physcons   import KCALMOL
from   Spline     import VadiSpline
#=======================================================#

NUMERR  = 0.1 # in kcal/mol

#=======================================================#
def ccvsic_checknfreqs(lcc_frqs,lic_frqs):
    if len(lic_frqs) == 0: return True
    idx = 0
    for ccfrqs,icfrqs in zip(lcc_frqs,lic_frqs):
        if len(ccfrqs) != len(icfrqs):
            return False
        idx += 1
    return True
#-------------------------------------------------------#
def ccvsic_checkts(data_x,lcc_frqs,lic_frqs):
    if len(lic_frqs) == 0: return True
    for idx,s_i in enumerate(data_x):
        if s_i != 0.0: continue
        ccfrqs = lcc_frqs[idx]
        icfrqs = lic_frqs[idx]
        if len(ccfrqs) != len(icfrqs): return False
        diff = max([abs(f1-f2) for f1,f2 in zip(ccfrqs,icfrqs)])*KCALMOL
        if diff > NUMERR: return False
        return True
#-------------------------------------------------------#
def ccvsic_checkzpe(lcc_tzpe,lic_tzpe):
    if len(lic_tzpe) == 0: return True
    for zpe_cc,zpe_ic in zip(lcc_tzpe,lic_tzpe):
        diff = (zpe_ic-zpe_cc)*KCALMOL
        if diff < -NUMERR: return False
    return True
#=======================================================#


def path2vadi(tcommon,drst,Eref=None,ics=None,boolint=False,lowfq={},freqscal=1.0):
    '''
    lowfq in case of imaginary frequencies
    '''
    if   boolint in [False,"no","No","n","N",None]: boolint = False
    elif boolint in [True,"yes","YES","y","Y"]    : boolint = True

    if boolint and ics in [None,False,[]]: raise Exc.NoICS(Exception)

    ch,mtp,atnums,masses,mu = tcommon

    # Sorted labels (by s)
    slabels = sd.sorted_points(drst,hess=True)

    # Reference energy
    if Eref is None:
       lbw, lfw, sbw, sfw, Ebw, Efw = sd.rstlimits(drst)
       Eref = Ebw

    # Independent variable
    data_x = [drst[label][0] for label in slabels]

    # mep energy
    listV0 = [drst[label][1]-Eref for label in slabels]

    # Dependent variable (cc)
    lcc_tzpe, lcc_frqs, lcc_Vadi = [], [], []
    lic_tzpe, lic_frqs, lic_Vadi = [], [], []
    dMols = {}
    for label in slabels:
        # data in drst
        s_i, E_i, xms_i,gms_i,Fms_i,v0_i,v1_i,t_i = drst[label]
        # project gradient
        if s_i == 0.0: bool_pg = False
        else         : bool_pg = True
        # lowfq
        if   s_i == 0: dlowfq = {}
        elif s_i  < 0: dlowfq = lowfq.get("bw",{})
        elif s_i  > 0: dlowfq = lowfq.get("fw",{})
        # mass-scaled --> Cartesian coords
        xcc = fncs.ms2cc_x(xms_i,masses,mu)
        gcc = fncs.ms2cc_g(gms_i,masses,mu)
        Fcc = fncs.ms2cc_F(Fms_i,masses,mu)
        # create Molecule instance
        mol = Molecule()
        mol.set(xcc,atnums,ch,mtp,E_i,gcc,Fcc,masses)
        mol.prep()
        mol.set_fscal(freqscal)
        mol.setup(mu,projgrad=bool_pg)
        mol.clean_freqs("cc")      # it may be needed with orca
        mol.deal_lowfq(dlowfq,"cc") # deal with low frequencies
        mol.ana_freqs("cc")        # calculate zpe
        # append data
        lcc_tzpe.append(float(mol._cczpe)   )
        lcc_Vadi.append(mol._ccV1 - Eref    )
        lcc_frqs.append(list(mol._ccfreqs))
        # internal coordinates
        if boolint:
           mol.icfreqs(ics,bool_pg)
           mol.deal_lowfq(dlowfq,"ic")
           mol.ana_freqs("ic")
           # append data
           lic_tzpe.append(float(mol._iczpe) )
           lic_Vadi.append(mol._icV1 - Eref  )
           lic_frqs.append(list(mol._icfreqs))
        # save instance
        dMols[label] = (s_i,mol)

    tuple_cc = (data_x,lcc_frqs,lcc_tzpe)
    tuple_ic = (data_x,lic_frqs,lic_tzpe)
    # Generate splines and setup
    Vadi_cc = VadiSpline(data_x,lcc_Vadi)
    if boolint: Vadi_ic = VadiSpline(data_x,lic_Vadi)
    else      : Vadi_ic = None
    # Return data
    return dMols, Vadi_cc, Vadi_ic, tuple_cc, tuple_ic, listV0
