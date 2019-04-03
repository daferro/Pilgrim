#!/usr/bin/python2.7
'''
*----------------------------------*
| Package    :  common             |
| Module     :  ispe               |
| Last Update:  2018/10/06 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

''' 

#==============================================#
import fncs      as     fncs
import physcons  as     pc
import numpy     as     np
import files     as     ff
import steepdesc as     sd
import Exceptions as    Exc
from   Molecule  import Molecule
from   Spline    import Spline
from   steepdesc import TSLABEL
from   criteria  import EPS_DLEVELS
#==============================================#

def read_ispe(filename):
    '''
    read input file for ispe
    '''
    # read lines
    lines = ff.read_file(filename)
    lines = fncs.clean_lines(lines,strip=True)
    # initialize data
    ispe_xy = []
    VR , VP = None, None
    tension = 0.0
    # find data in lines
    for line in lines:
        if line == "\n": continue
        label, val = line.split()
        val = float(val)
        if   label.lower() == "tension": tension = val
        elif label.lower() == "reac"   : VR = val
        elif label.lower() == "prod"   : VP = val
        else                           : ispe_xy.append( (label,val) )
    return ispe_xy, tension, VR, VP

#==============================================#
def get_s0_L(VR,VTS,VP,mu,ifreq):
    '''
    calculates s0 and L in base of
    low-level parameters
    '''
    sA0 = -np.sqrt( (VTS-VR)/(ifreq**2)/mu )
    sB0 = +np.sqrt( (VTS-VP)/(ifreq**2)/mu )
    
    sA  = -min(  abs(sA0),2*sB0)
    sB  = +min(2*abs(sA0),  sB0)

    s0 = (sA+sB)/2.0
    L  = (abs(sA)+sB)/2.0
    return s0, L
#----------------------------------------------#
def s2z(sval,s0,L):
    '''
    converts s --> z 
    '''
    arg = (sval-s0)/L
    z  = 2.0/pc.PI * np.arctan(arg)
    return z
#----------------------------------------------#
def get_ts_ifreq(xcc,gcc,Fcc,E,tcommon):
    ch, mtp, atonums, masses, mu = tcommon
    ts = Molecule()
    ts.set(xcc,atonums,ch,mtp,E,gcc,Fcc,masses)
    ts.prep()
    ts.setup(mu)
    ts.ana_freqs()
    [(ifreq, evec)] = ts._ccimag
    return ifreq, mu
#----------------------------------------------#
def in_interval(si,s0,sn):
    return s0-EPS_DLEVELS <= si <= sn+EPS_DLEVELS
#----------------------------------------------#
def gen_data_for_spline(ispe_xy,drst,s0,L):
    lz, ldE = [], []
    for sval,point,Vhl in ispe_xy:
        # low-level data
        Vll = drst[point][1]
        # convert to z-val
        z = s2z(sval,s0,L)
        # get difference
        dE = Vll-Vhl
        # append data
        lz.append(z)
        ldE.append(dE)
    return lz, ldE
#==============================================#

def ispe(tcommon,drst,ispe_xy,tension):
    '''
    V^LL --> V^HL
    '''

    # points in rst
    points = sd.sorted_points(drst,hess=False)

    # Get imaginary frequency for TS (low-level)
    sTS, ETS_ll, xcc, gcc, Fcc = drst[sd.TSLABEL][0:5]
    ifreq, mu = get_ts_ifreq(xcc,gcc,Fcc,ETS_ll,tcommon)
    
    # convert x in ispe_xy to labels
    for idx, (l_i,V_i) in enumerate(ispe_xy):
        if l_i in points:
           s_i = drst[l_i][0]
           label_i = l_i
        else:
           s_i = float(l_i)
           label_i = None
           for point in points:
               if abs(drst[point][0]-s_i) < EPS_DLEVELS: label_i = point
           if label_i is None: raise Exc.DLEVELsthWrong(Exception)
        ispe_xy[idx] = (s_i,label_i,V_i)
    # sort points
    ispe_xy.sort()

    # only one point
    if len(ispe_xy) == 1:
       sx, lx, Ex_hl = ispe_xy[0]
       Ex_ll = drst[lx][1]
       # calculate energy difference
       diffE = Ex_ll - Ex_hl
       # points
       points = sd.sorted_points(drst,hess=False)
       # save old energies
       xx   = [ drst[point][0] for point in points]
       yyll = [ drst[point][1] for point in points]
       # apply difference to all points
       for point in points:
           drst[point] = list(drst[point])
           drst[point][1] = drst[point][1] - diffE
           drst[point] = tuple(drst[point])
    # more than one point
    else:
       # reduce points of drst to those for DLEVEL
       s1, l1, E1_hl = ispe_xy[0]
       sn, ln, En_hl = ispe_xy[-1]
       drst   = {point:drst[point] for point in points if in_interval(drst[point][0],s1,sn)}
       points = sd.sorted_points(drst,hess=False)

       # get s0 and L from lower level
       E1_ll, En_ll = drst[l1][1], drst[ln][1]
       s0, L = get_s0_L(E1_ll,ETS_ll,En_ll,mu,ifreq)

       # generate data for spline
       lxx, lyy = gen_data_for_spline(ispe_xy,drst,s0,L)

       # create spline
       spl = Spline(lxx,lyy,tension=tension)

       # save old energies
       xx   = [ drst[point][0] for point in points]
       yyll = [ drst[point][1] for point in points]

       # modify drst
       for point in points:
           s, Vll =  drst[point][0:2]
           z      = s2z(s,s0,L)
           drst[point] = list(drst[point])
           drst[point][1] = Vll - spl(z)
           drst[point] = tuple(drst[point])

    # save new energies
    yyhl = [ drst[point][1] for point in points]

    # return data
    return drst, points, xx, yyll, yyhl

#def ispe(tcommon,drst,VR_ll,VP_ll,ispe_inp):
#    '''
#    drst   : dict with the MEP
#    list_lb: labels of points with HL calculations
#    list_HL: energies at HL
#    '''
#    # read data in input
#    ispe_xy, tension, VR_hl, VP_hl = read_ispe(ispe_inp)
#
#    # Get imaginary frequency for TS (low-level)
#    sTS, VTS_ll, xcc, gcc, Fcc = drst[sd.TSLABEL][0:5]
#    ifreq, mu = get_ts_ifreq(xcc,gcc,Fcc,VTS_ll,tcommon)
#    
#    # Energy difference in reactants
#    ER_diff = VR_ll-VR_hl
#
#    # get s0 and L from lower level
#    s0, L = get_s0_L(VR_ll,VTS_ll,VP_ll,mu,ifreq)
#
#    # generate data for spline
#    lxx, lyy = gen_data_for_spline(ispe_xy,drst,s0,L)
#
#    # add data for reactant and product
#    dE = VR_hl-VR_ll
#    lxx = [-1.0]+lxx
#    lyy = [ dE ]+lyy
#    dE = VP_hl-VP_ll
#    lxx = lxx+[+1.0]
#    lyy = lyy+[ dE ]
#
#    # correct energy difference (same origin LL and HL)
#    lyy = [Ei+ER_diff for Ei in lyy]
#
#    # create spline
#    spl = Spline(lxx,lyy,tension=tension)
#
#    # save old energies
#    points = sd.sorted_points(drst,hess=False)
#    xx   = [ drst[point][0] for point in points]
#    yyll = [ drst[point][1] for point in points]
#
#    # correct energies
#    for point in points:
#        s, Vll =  drst[point][0:2]
#        z      = s2z(s,s0,L)
#        drst[point] = list(drst[point])
#        drst[point][1] = Vll + spl(z)
#        drst[point] = tuple(drst[point])
#
#    # save new energies
#    yyhl = [ drst[point][1] for point in points]
#
#    # return data
#    return drst, xx, yyll, yyhl
#
#
#def test():
#    VR_ll = -186.21788
#    VP_ll = -186.21081
#    rstfile  = "/home/david/PyFerro/hcooh/IOfiles/RST/TS.001.rst"
#    ispe_inp = "/home/david/PyFerro/hcooh/ispe.inp"
#    tpath, tcommon, drst = ff.read_rst(rstfile)
#    drst, xx, yyll, yyhl = ispe(tcommon,drst,VR_ll,VP_ll,ispe_inp)
#    import matplotlib.pyplot as plt
#    plt.plot(xx,yyll,'k--')
#    plt.plot(xx,yyhl,'r--')
#    plt.show()

#if __name__ == "__main__": test()
