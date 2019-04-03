#!/usr/bin/python2.7
'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 1.0
License     : MIT/x11

Copyright (c) 2019, David Ferro Costas (david.ferro@usc.es) and
Antonio Fernandez Ramos (qf.ramos@usc.es)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
---------------------------

*----------------------------------*
| Package    :  common             |
| Module     :  internal           |
| Last Update:  2019/04/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains some functions related
to internal coordinates

Functions (internal.py):
  * ic2string(ic)
  * string2ic(icstring)
  * merge_ics(ics_st,ics_ab,ics_lb,ics_it,ics_pt)
  * unmerge_ics(all_ics)
  * count_ics(all_ics)
  * get_adjmatrix(xcc,symbols,scale=1.2,mode="bool")
  * get_numbonds(amatrix)
  * adjacency_matrix2list(amatrix)
  * get_subgraph(node,alist,fragment=[])
  * get_fragments(alist)
  * distance_2fragments(frg1,frg2,xcc)
  * distance_allfragments(fragments,xcc)
  * frags_distances(fragments)
  * link_fragments(xcc,amatrix,nfrags=1)
  * ics_value(xcc,ic)
  * ics_get_stretchings(cmatrix,natoms)
  * ics_get_iccentral(adj_list)
  * ics_classify_bends(xcc,ic_3ats,eps_lin=5.0)
  * ics_get_ptorsions(ics_st,adj_list,xcc,epslin=5.0)
  * ics_get_ltorsions(ics_lb,adj_list)
  * ics_from_geom(xcc,symbols,scale=1.3,nfrags=1,eps_lin=5.0)
  * ics_from_gts(gtsfile,scale=1.3,nfrags=1,eps_lin=5.0)
  * ics_depure_bendings(ic_3ats,keep=[])
  * ics_depure_itorsions(ic_4ats,keep=[])
  * ics_depure_ptorsions(ic_4ats,keep=[])
  * ics_depure(ics,keep=[])
  * ics_correctdir(x1,evec,ic,sign)
  * ics_idir(xcc,symbols,masses,freqs,ms_evecs,ics=[],mu=1.0/AMU)
  * wilson_bvecs(xcc)
  * wilson_stretch(bond_vectors,ij,natoms)
  * wilson_abend(bond_vectors,ijk,natoms)
  * wilson_auxlinB(m,o,n,k)
  * wilson_auxlinC(m,o,n,k,Bk=None,Dk=None)
  * wilson_lbend(m,o,n,MON,natoms)
  * wilson_torsion(bond_vectors,ijkl,natoms)
  * wilson_getBC(xcc,all_ics)
  * wilson_getu(masses)
  * wilson_getG(u,B)
  * wilson_gf_internal(u,B,C,Ginv,gcc,Fcc)
  * wilson_gf_nonred(G,Ginv,g,f)
  * wilson_prj_rc(gnr,fnr,G,nics)
  * wilson_evecsincart(L,G,A,masses)
  * calc_icfreqs(Fcc,masses,xcc,gcc,all_ics,bool_prc=False)
  * nonredundant(xcc,masses,gcc,Fcc,all_ics,ccfreqs,unremov=[],ncycles=None)
  * nonredundant_gtsfiles(gtsfiles,all_ics,unremov=[],ncycles=None)
'''

import random
import os
import sys
import numpy          as     np
import fncs           as     fncs
from   scipy.optimize import brenth
from   dicts          import dpt_s2cr
from   physcons       import AMU
from   files          import read_gtsfile
from   criteria       import EPS_SCX
from   criteria       import EPS_SVD
from   criteria       import EPS_GIV
from   criteria       import EPS_ICF
from   criteria       import EPS_SMALLANGLE

#===============================================#
# Internal coordinates / Graph Theory           #
#===============================================#
def ic2string(ic):
    ictype, icatoms = ic
    if ictype == "st": return "-".join(["%i"%(at+1) for at in icatoms])
    if ictype == "ab": return "-".join(["%i"%(at+1) for at in icatoms])
    if ictype == "pt": return "-".join(["%i"%(at+1) for at in icatoms])
    if ictype == "lb": return "=".join(["%i"%(at+1) for at in icatoms])
    if ictype == "it": return "_".join(["%i"%(at+1) for at in icatoms])
#---------------------------------------------#
def string2ic(icstring):
    if "-" in icstring:
        atoms = [int(at)-1 for at in icstring.split("-")]
        if   len(atoms) == 2: case = "st"
        elif len(atoms) == 3: case = "ab"
        elif len(atoms) == 4: case = "pt"
        else: exit("Problems with internal coordinate!")
        if atoms[0] > atoms[-1]: atoms = atoms[::-1]
    if "=" in icstring:
        atoms = [int(at)-1 for at in icstring.split("=")]
        case  = "lb"
        if len(atoms) != 3: exit("Problems with internal coordinate!")
        if atoms[0] > atoms[-1]: atoms = atoms[::-1]
    if "_" in icstring:
        atoms = [int(at)-1 for at in icstring.split("_")]
        case  = "it"
        if len(atoms) != 4: exit("Problems with internal coordinate!")
        atoms = tuple(sorted(atoms[0:3])+atoms[3:4])
    return (case,atoms)
#-----------------------------------------------#
def merge_ics(ics_st,ics_ab,ics_lb,ics_it,ics_pt):
    all_ics  = [ ("st",ic) for ic in sorted(ics_st)] # stretching
    all_ics += [ ("ab",ic) for ic in sorted(ics_ab)] # angular bending
    all_ics += [ ("lb",ic) for ic in sorted(ics_lb)] # linear  bending
    all_ics += [ ("it",ic) for ic in sorted(ics_it)] # improper torsion
    all_ics += [ ("pt",ic) for ic in sorted(ics_pt)] # proper   torsion
    return all_ics
#-----------------------------------------------#
def unmerge_ics(all_ics):
    ics_st = [ ic for ic_type,ic in all_ics if ic_type=="st"]
    ics_ab = [ ic for ic_type,ic in all_ics if ic_type=="ab"]
    ics_lb = [ ic for ic_type,ic in all_ics if ic_type=="lb"]
    ics_it = [ ic for ic_type,ic in all_ics if ic_type=="it"]
    ics_pt = [ ic for ic_type,ic in all_ics if ic_type=="pt"]
    return ics_st,ics_ab,ics_lb,ics_it,ics_pt
#-----------------------------------------------#
def count_ics(all_ics):
    ics_st,ics_ab,ics_lb,ics_it,ics_pt = unmerge_ics(all_ics)
    nICs = len(ics_st)+len(ics_ab)+2*len(ics_lb)+len(ics_it)+len(ics_pt)
    return nICs
#-----------------------------------------------#
def get_adjmatrix(xcc,symbols,scale=1.2,mode="bool"):
    '''
    returns adjacency matrix (connection matrix);
    also distance matrix and number of bonds
    * mode = bool, int
    '''
    nbonds  = 0
    nat     = fncs.howmanyatoms(xcc)
    dmatrix = fncs.get_distmatrix(xcc)
    if mode == "bool": no, yes = False, True
    if mode == "int" : no, yes = 0    , 1
    cmatrix = [ [no for ii in range(nat)] for jj in range(nat)]
    for ii in range(nat):
        cr_ii = dpt_s2cr[symbols[ii]] # covalent radius
        for jj in range(ii+1,nat):
            cr_jj = dpt_s2cr[symbols[jj]]
            dref  = (cr_ii+cr_jj)*scale
            if dmatrix[ii][jj] < dref:
               nbonds += 1
               cmatrix[ii][jj] = yes
               cmatrix[jj][ii] = yes
    return cmatrix, dmatrix, nbonds
#-----------------------------------------------#
def get_numbonds(amatrix):
    nbonds = 0
    nnodes = len(amatrix)
    for node1 in range(nnodes):
        for node2 in range(node1+1,nnodes):
            if amatrix[node1][node2] in [True,1]: nbonds += 1
    return nbonds
#-----------------------------------------------#
def adjacency_matrix2list(amatrix):
    alist = {}
    for node1,row in enumerate(amatrix):
       alist[node1]= [node2 for node2, bonded in enumerate(row) if bonded in [True,1]]
    return alist
#-----------------------------------------------#
def get_subgraph(node,alist,fragment=[]):
    fragment += [node]
    neighbors = alist[node]
    for neighbor in neighbors:
        if neighbor not in fragment:
            fragment = get_subgraph(neighbor,alist,fragment)
    return fragment
#-----------------------------------------------#
def get_fragments(alist):
    fragments = []
    visited   = set([])
    for node in alist.keys():
        if node in visited: continue
        fragment = set(get_subgraph(node,alist,[]))
        if fragment not in fragments: fragments.append(fragment)
        visited = visited.union(fragment)
    return fragments
#-----------------------------------------------#
def distance_2fragments(frg1,frg2,xcc):
    min_dist = float("inf")
    pair     = (None,None)
    for at1 in frg1:
        x1 = fncs.xyz(xcc,at1)
        for at2 in frg2:
            x2 = fncs.xyz(xcc,at2)
            dist = fncs.distance(x1,x2)
            if dist < min_dist:
               min_dist = dist
               pair     = (at1,at2)
    return min_dist, pair
#-----------------------------------------------#
def distance_allfragments(fragments,xcc):
    nfrags = len(fragments)
    the_list = []
    for idx1 in range(nfrags):
        frg1 = fragments[idx1]
        for idx2 in range(idx1+1,nfrags):
            frg2 = fragments[idx2]
            dist, (at1,at2) = distance_2fragments(frg1,frg2,xcc)
            the_list.append( (dist,idx1,idx2,at1,at2) )
    return sorted(the_list)
#-----------------------------------------------#
def frags_distances(fragments):
    ''' use distance_allfragments -  this one returns distance
        of 1 and just first atom in each fragment'''
    fdists = []
    nfrags = len(fragments)
    for idx1 in range(nfrags):
        for idx2 in range(idx1+1,nfrags):
            dist = 1.0
            atf1 = list(fragments[idx1])[0]
            atf2 = list(fragments[idx2])[0]
            fdists.append( (1.0,idx1,idx2,atf1,atf2) )
    fdists.sort()
    return fdists
#-----------------------------------------------#
def link_fragments(xcc,amatrix,nfrags=1):
    if   amatrix[0][0] is False: bonded = True
    elif amatrix[0][0] is 0    : bonded = 1
    else: exit("sth wrong in adjacency matrix!")
    alist     = adjacency_matrix2list(amatrix)
    fragments = get_fragments(alist)
   #fdists    = frags_distances(fragments)
    fdists    = distance_allfragments(fragments,xcc)
    inumfrags = len(fragments)
    fnumfrags = len(fragments)
    for dist,idx1,idx2,atf1,atf2 in fdists:
        fragments[idx1] = fragments[idx1].union(fragments[idx2])
        fragments[idx2] = set([])
        amatrix[atf1][atf2] = bonded
        amatrix[atf2][atf1] = bonded
        fnumfrags = sum([1 for frag in fragments if len(frag)!=0])
        if fnumfrags == nfrags: break
    fragments = [frag for frag in fragments if len(frag) != 0]
    return amatrix, fragments, inumfrags, fnumfrags
#-----------------------------------------------#
def ics_value(xcc,ic):
    ''' ic = (ic_type,ic_atoms)'''
    if type(ic) == type("string"):
        ic_type,ic_atoms = string2ic(ic)
    else:
        ic_type,ic_atoms = ic
    if ic_type == "st": return fncs.distance( *(fncs.xyz(xcc,at) for at in ic_atoms) )
    if ic_type == "ab": return fncs.angle(    *(fncs.xyz(xcc,at) for at in ic_atoms) )
    if ic_type == "lb": return fncs.angle(    *(fncs.xyz(xcc,at) for at in ic_atoms) )
    if ic_type == "it": return fncs.dihedral( *(fncs.xyz(xcc,at) for at in ic_atoms) )
    if ic_type == "pt": return fncs.dihedral( *(fncs.xyz(xcc,at) for at in ic_atoms) )
#-----------------------------------------------#
def ics_get_stretchings(cmatrix,natoms):
   ics_st = [(at1,at2) for at1 in range(natoms) for at2 in range(at1,natoms) if cmatrix[at1][at2] in [True,1]]
   return ics_st
#-----------------------------------------------#
def ics_get_iccentral(adj_list):
    ic_3ats = []
    ic_4ats = []
    for at2 in adj_list.keys():
        bonded  = adj_list[at2]
        nbonded = len(bonded)
        if nbonded < 2: continue
        for idx1 in range(nbonded):
            for idx3 in range(idx1+1,nbonded):
                bending = (bonded[idx1],at2,bonded[idx3])
                ic_3ats.append(bending)
                if nbonded < 3: continue
                for idx4 in range(idx3+1,nbonded):
                    improper_torsion = (bonded[idx1],bonded[idx3],bonded[idx4],at2)
                    ic_4ats.append(improper_torsion)
    return ic_3ats, ic_4ats
#-----------------------------------------------#
def ics_classify_bends(xcc,ic_3ats,eps_lin=5.0):
    ics_lb = []
    ics_ab = []
    thetas     = {}
    for at1,at2,at3 in ic_3ats:
       x1 = fncs.xyz(xcc,at1)
       x2 = fncs.xyz(xcc,at2)
       x3 = fncs.xyz(xcc,at3)
       theta = abs(fncs.rad2deg(fncs.angle(x1,x2,x3)))
       if theta < eps_lin or theta > 180-eps_lin:
          ics_lb.append( (at1,at2,at3) )
       else:
          ics_ab.append( (at1,at2,at3) )
       thetas[(at1,at2,at3)] = theta
    return ics_lb, ics_ab, thetas
#-----------------------------------------------#
def ics_get_ptorsions(ics_st,adj_list,xcc,epslin=5.0):
    '''epslin in degrees'''
    ics_pt = []
    for at2,at3 in ics_st:
        x2 = fncs.xyz(xcc,at2)
        x3 = fncs.xyz(xcc,at3)
        bondedto2 = list(adj_list[at2])
        bondedto3 = list(adj_list[at3])
        bondedto2.remove(at3)
        bondedto3.remove(at2)
        if len(bondedto2) == 0: continue
        if len(bondedto3) == 0: continue
        for at1 in bondedto2:
            x1 = fncs.xyz(xcc,at1)
            for at4 in bondedto3:
                if at1 == at4: continue
                x4 = fncs.xyz(xcc,at4)
                # the two angles
                angA  = fncs.rad2deg(fncs.angle(x1,x2,x3))
                angB  = fncs.rad2deg(fncs.angle(x2,x3,x4))
                # linear?
                booleanA = angA < epslin or angA > 180-epslin
                booleanB = angB < epslin or angB > 180-epslin
                if booleanA or booleanB: continue
                ptorsion = (at1,at2,at3,at4)
                ics_pt.append(ptorsion)
    return ics_pt
#-----------------------------------------------#
def ics_get_ltorsions(ics_lb,adj_list):
    ic_ltors = []
    for at1,at2,at3 in ics_lb:
        bondedto1 = list(adj_list[at1])
        bondedto3 = list(adj_list[at3])
        if at1 in bondedto3: bondedto3.remove(at1)
        if at2 in bondedto3: bondedto3.remove(at2)
        if at2 in bondedto1: bondedto1.remove(at2)
        if at3 in bondedto1: bondedto1.remove(at3)
        for at0 in bondedto1:
            for at4 in bondedto3:
                if at0 == at4: continue
                ltorsion = (at0,at1,at3,at4)
                ic_ltors.append(ltorsion)
    return ic_ltors
#-----------------------------------------------#
def ics_from_geom(xcc,symbols,scale=1.3,nfrags=1,eps_lin=5.0):
    natoms = len(symbols)
    amatrix, dmatrix, nbonds = get_adjmatrix(xcc,symbols,scale=scale,mode="bool")
    amatrix, fragments, inumfrags, fnumfrags = link_fragments(xcc,amatrix,nfrags=nfrags)
    if inumfrags != fnumfrags:
       nbonds = get_numbonds(amatrix)
    alist = adjacency_matrix2list(amatrix)
    ics_st = ics_get_stretchings(amatrix,natoms)
    ic_3ats, ics_it = ics_get_iccentral(alist)
    ics_lb, ics_ab, angles = ics_classify_bends(xcc,ic_3ats,eps_lin)
    ics_pt  = ics_get_ptorsions(ics_st,alist,xcc,eps_lin)
    ics_pt += ics_get_ltorsions(ics_lb,alist)
    return merge_ics(ics_st,ics_ab,ics_lb,ics_it,ics_pt)
#-----------------------------------------------#
def ics_from_gts(gtsfile,scale=1.3,nfrags=1,eps_lin=5.0):
    xcc,atonums,ch,mtp,E,gcc,Fcc,masses,pgroup,rotsigma,freq_list = read_gtsfile(gtsfile)
    symbols = fncs.get_symbols(atonums)
    return ics_from_geom(xcc,symbols,scale,nfrags,eps_lin)
#-----------------------------------------------#
def ics_depure_bendings(ic_3ats,keep=[]):
    ''' up to 3 bendings per central atom'''
    centers = [at2 for at1,at2,at3 in keep]
    random.shuffle(ic_3ats)
    for at1,at2,at3 in ic_3ats:
        if centers.count(at2) >= 3: continue
        #if at2 in centers: continue
        keep.append( (at1,at2,at3) )
        centers.append(at2)
    return keep
#-----------------------------------------------#
def ics_depure_itorsions(ic_4ats,keep=[]):
    ''' up to two improper torsions per atom'''
    centers = [at4 for at1,at2,at3,at4 in keep]
    random.shuffle(ic_4ats)
    for at1,at2,at3,at4 in ic_4ats:
        if centers.count(at4) >= 2: continue
        #if at4 in centers: continue
        keep.append( (at1,at2,at3,at4) )
        centers.append(at4)
    return keep
#-----------------------------------------------#
def ics_depure_ptorsions(ic_4ats,keep=[]):
    ''' one torsion per bond'''
    centers = [(at2,at3) for at1,at2,at3,at4 in keep]
    random.shuffle(ic_4ats)
    for at1,at2,at3,at4 in ic_4ats:
        if (at2,at3) in centers: continue
        if (at3,at2) in centers: continue
        keep.append( (at1,at2,at3,at4) )
        centers.append((at2,at3))
    return keep
#-----------------------------------------------#
def ics_depure(ics,keep=[]):
    ''' keep: those that cannot be removed '''
    ics_st ,ics_ab ,ics_lb ,ics_it ,ics_pt  = unmerge_ics(ics)
    keep_st,keep_ab,keep_lb,keep_it,keep_pt = unmerge_ics(keep)
    # depure angular bendings
    ics_ab = ics_depure_bendings(ics_ab,keep_ab)
    # depure torsions
    ics_it = ics_depure_itorsions(ics_it,keep_it)
    ics_pt = ics_depure_ptorsions(ics_pt,keep_pt)
    # merge again
    ics = merge_ics(ics_st,ics_ab,ics_lb,ics_it,ics_pt)
    return ics
#-----------------------------------------------#
def ics_correctdir(x1,evec,ic,sign,masses=None,mu=None):
    '''
    x1 and evec NOT in mass-scaled
    '''
    x2 = [xi+ei for xi,ei in zip(x1,evec)]
    if masses is not None:
       x1 = fncs.ms2cc_x(x1,masses,mu)
       x2 = fncs.ms2cc_x(x2,masses,mu)
    val1 = ics_value(x1,ic)
    val2 = ics_value(x2,ic)
    diff = val2-val1
    if   diff > 0.0:
        if sign == "++": return True
        if sign == "--": return False
    elif diff < 0.0:
        if sign == "++": return False
        if sign == "--": return True
#-----------------------------------------------#
def ics_idir(xcc,symbols,masses,freqs,ms_evecs,ics=[],mu=1.0/AMU):
    '''
    returns the IC which varies the most due to the imaginary frequency
    '''
    if len(ics) == 0: ics = ics_from_geom(xcc,symbols)
    for freq,Lms in zip(freqs,ms_evecs):
        # only imaginary frequency
        if freq >= 0.0: continue
        Lcc = fncs.ms2cc_x(Lms,masses,mu)
        xfin = [xi+ei for xi,ei in zip(xcc,Lcc)]
        target_ic   = None
        target_sign = None
        maxdiff     = -float("inf")
        for ic in ics:
            ival = ics_value(xcc ,ic)
            fval = ics_value(xfin,ic)
            # the sign
            if fval >= ival: sign = "++"
            else           : sign = "--"
            # reference for bonds or angles
            if len(ic[1]) == 2: reference = 1.0       # 1.0 bohr
            else              : reference = np.pi/2.0 # 90 degrees
            # get absolute diff
            if len(ic[1]) == 4: adiff = abs(fncs.angular_dist(fval,ival,'rad'))
            else              : adiff = abs(fval - ival)
            # get relative difference with regards to reference
            reldiff = abs(adiff/reference)
            # get the one that changes the most
            if reldiff > maxdiff:
                target_ic   = ic
                target_sign = sign
                maxdiff     = reldiff
        return target_ic,target_sign
#===============================================#


#===============================================#
# WILSON method for freqs in internal coords    #
#===============================================#
'''
Some references:
 [1] D.F. McIntosh and K.H. Michaelian, The Wilson GF Matrix Method of Vibrational Analysis. Part I, II, and III
     Can.J.Spectros., 24, 1-10 (1979) ; Can.J.Spectros., 24, 35-40 (1979) ; Can.J.Spectros., 24, 65-74 (1979)
 [2] Ian H. Williams, Torsional Internal Coordinates in Normal Coordinate Calculations
     J. Mol. Spectros., 66, 288-301 (1977)
 [3] C.F. Jackels, Z. Gu, D.G. Truhlar, Reaction-path potential and vibrational frequencies in terms of curvilinear internal coordinates
     J. Chem. Phys. 102, 3188-3201 (1995)
 [4] Y-Y Chuang and D.G. Truhlar, Reaction-Path Dynamics in Redundant Internal Coordinates
     J. Phys. Chem. A, 102, 242-247 (1998)
 [5] V. Bakken and T. Helgaker, The efficient optimization of molecular geometries using redundant internal coordinates
     J. Chem. Phys., 117, 9160-9174 (2002)
'''
#-----------------------------------------------#
def wilson_bvecs(xcc):
    '''
    This functions calculates the distance
    between each pair of atoms (dij) and also
    the unit bond vector eij = (rj - ri)/dij
    '''
    nat = fncs.howmanyatoms(xcc)
    bond_vectors = {}
    for ii in range(nat):
        ri = np.array(fncs.xyz(xcc,ii))
        for jj in range(ii+1,nat):
            rj = np.array(fncs.xyz(xcc,jj))
            eij = rj-ri
            dij = np.linalg.norm(eij)
            eij = eij / dij
            bond_vectors[(ii,jj)] = ( eij, dij)
            bond_vectors[(jj,ii)] = (-eij, dij)
    return bond_vectors
#-----------------------------------------------#
def wilson_stretch(bond_vectors,ij,natoms):
    '''
    Returns the row of the B matrix associated to
    the bond length between atom i and atom j and
    also the corresponding C matrix.
    Check ref [1] and [5].
    '''
    # Using nomenclature of reference [5]
    n = min(ij)
    m = max(ij)
    u, r = bond_vectors[(n,m)]
    # Calculating 1st derivatives: row of B matrix
    B_row = [0.0 for idx in range(3*natoms)]
    for a in [m,n]:
        if a == m: zeta_amn = +1.0
        if a == n: zeta_amn = -1.0
        for i in [0,1,2]:
            dr_dai = zeta_amn * u[i]
            B_row[3*a+i] = dr_dai
    # Calculating 2nd derivatives: 2D matrix of C tensor
    C_matrix = [ [0.0 for idx1 in range(3*natoms)] for idx2 in range(3*natoms)]
    #C_matrix = np.zeros( (3*natoms,3*natoms) )
    for a in [m,n]:
      for i in [0,1,2]:
          for b in [m,n]:
            for j in [0,1,2]:
               if C_matrix[3*a+i][3*b+j] != 0.0: continue
               # Get delta values
               if a == b: delta_ab = 1.0
               else:      delta_ab = 0.0
               if i == j: delta_ij = 1.0
               else:      delta_ij = 0.0
               # Get C element
               dr_daidbj = ((-1.0) ** delta_ab) * (u[i]*u[j] - delta_ij) / r
               # Save data in both positions
               C_matrix[3*a+i][3*b+j] = dr_daidbj
               C_matrix[3*b+j][3*a+i] = dr_daidbj
    return [B_row], [C_matrix]
#-----------------------------------------------#
def wilson_abend(bond_vectors,ijk,natoms):
    '''
    Returns the row of the B matrix associated to the
    i-j-k angle bend and the corresponding C matrix.
    Check ref [1] and [5].
    '''
    # Using nomenclature of reference [5]
    m, o, n = ijk
    u, lambda_u = bond_vectors[(o,m)]
    v, lambda_v = bond_vectors[(o,n)]
    # Get internal coordinate: bond angle
    q = fncs.angle_vecs(u,v)
    sinq = np.sin(q)
    cosq = np.cos(q)
    # Generation of w
    w = np.cross(u,v)
    w = w / np.linalg.norm(w)
    uxw = np.cross(u,w)
    wxv = np.cross(w,v)
    # Calculating 1st derivatives: row of B matrix
    B_row = [0.0 for idx in range(3*natoms)]
    for a in [m,o,n]:
        # Get zeta values
        if a == m: zeta_amo = +1.0; zeta_ano =  0.0
        if a == o: zeta_amo = -1.0; zeta_ano = -1.0
        if a == n: zeta_amo =  0.0; zeta_ano = +1.0
        for i in [0,1,2]:
            # Get B element
            dq_dai =  zeta_amo * uxw[i] / lambda_u + zeta_ano * wxv[i] / lambda_v
            B_row[3*a+i] = dq_dai
    # Calculating 2nd derivatives: 2D matrix of C tensor
    #C_matrix = np.zeros( (3*natoms,3*natoms) )
    C_matrix = [ [0.0 for idx1 in range(3*natoms)] for idx2 in range(3*natoms)]
    if abs(sinq) < EPS_SCX: return [B_row], [C_matrix]
    for a in [m,o,n]:
      for i in [0,1,2]:
        for b in [m,o,n]:
          for j in [0,1,2]:
            if C_matrix[3*a+i][3*b+j] != 0.0: continue
            # Define all delta and zeta values
            if a == m: zeta_amo = +1.0; zeta_ano =  0.0
            if a == o: zeta_amo = -1.0; zeta_ano = -1.0
            if a == n: zeta_amo =  0.0; zeta_ano = +1.0
            if b == m: zeta_bmo = +1.0; zeta_bno =  0.0
            if b == o: zeta_bmo = -1.0; zeta_bno = -1.0
            if b == n: zeta_bmo =  0.0; zeta_bno = +1.0
            if i == j: delta_ij = 1.0
            else:      delta_ij = 0.0
            # Get second derivative
            t1 = zeta_amo*zeta_bmo*(u[i]*v[j]+u[j]*v[i]-3*u[i]*u[j]*cosq+delta_ij*cosq)/(lambda_u**2 * sinq)
            t2 = zeta_ano*zeta_bno*(v[i]*u[j]+v[j]*u[i]-3*v[i]*v[j]*cosq+delta_ij*cosq)/(lambda_v**2 * sinq)
            t3 = zeta_amo*zeta_bno*(u[i]*u[j]+v[j]*v[i]-u[i]*v[j]*cosq-delta_ij)/(lambda_u*lambda_v*sinq)
            t4 = zeta_ano*zeta_bmo*(v[i]*v[j]+u[j]*u[i]-v[i]*u[j]*cosq-delta_ij)/(lambda_u*lambda_v*sinq)
            t5 = cosq / sinq * B_row[3*a+i] * B_row[3*b+j]
            dr_daidbj = t1 + t2 + t3 + t4 - t5
            C_matrix[3*a+i][3*b+j] = dr_daidbj
            C_matrix[3*b+j][3*a+i] = dr_daidbj
    return [B_row], [C_matrix]
#-----------------------------------------------#
def wilson_auxlinB(m,o,n,k):
        om = m - o
        on = n - o
        dom = np.linalg.norm(om)
        don = np.linalg.norm(on)
        qk = ((n[2]-o[2])*(m[k]-o[k]) - (n[k]-o[k])*(m[2]-o[2])) / dom / don
        Bk = []
        Dk = []
        for a in ["m","o","n"]:
            for i in [0,1,2]:
                if a == "m": dam, dao, dan = 1.0, 0.0, 0.0
                if a == "o": dam, dao, dan = 0.0, 1.0, 0.0
                if a == "n": dam, dao, dan = 0.0, 0.0, 1.0
                if i ==  2 : di2 = 1.0
                else       : di2 = 0.0
                if i == k  : dik = 1.0
                else       : dik = 0.0
                # Numerator
                N_ai = dik*(dao-dan)*m[2] + di2*(dan-dao)*m[k] +\
                       dik*(dan-dam)*o[2] + di2*(dam-dan)*o[k] +\
                       dik*(dam-dao)*n[2] + di2*(dao-dam)*n[k]
                # Denominator
                D_ai = (dam-dao)*(m[i]-o[i]) * don / dom + \
                       (dan-dao)*(n[i]-o[i]) * dom / don
                Dk.append(D_ai)
                # Whole derivative
                B_ai = (N_ai - qk*D_ai) / (dom*don)
                Bk.append(B_ai)
        return np.array(Bk), np.array(Dk)
#-----------------------------------------------#
def wilson_auxlinC(m,o,n,k,Bk=None,Dk=None):
        om = m - o
        on = n - o
        dom = np.linalg.norm(om)
        don = np.linalg.norm(on)
        if (Bk is None) or (Dk is None): Bk, Dk = get_B(m,o,n,k,True)
        Ck = np.zeros( (9,9) )
        for a in ["m","o","n"]:
            for i in [0,1,2]:
                if a == "m": dam, dao, dan = 1.0, 0.0, 0.0; A=0
                if a == "o": dam, dao, dan = 0.0, 1.0, 0.0; A=1
                if a == "n": dam, dao, dan = 0.0, 0.0, 1.0; A=2
                if i ==  2 : di2 = 1.0
                else       : di2 = 0.0
                if i == k  : dik = 1.0
                else       : dik = 0.0
                for b in ["m","o","n"]:
                    for j in [0,1,2]:
                        if b == "m": dbm, dbo, dbn = 1.0, 0.0, 0.0; B=0
                        if b == "o": dbm, dbo, dbn = 0.0, 1.0, 0.0; B=1
                        if b == "n": dbm, dbo, dbn = 0.0, 0.0, 1.0; B=2
                        if j ==  2 : dj2 = 1.0
                        else       : dj2 = 0.0
                        if j == k  : djk = 1.0
                        else       : djk = 0.0
                        # Term 1
                        term1 = (djk*di2-dik*dj2)*( (dan-dao)*dbm + (dam-dan)*dbo + (dao-dam)*dbn )
                        # Term 2
                        term2 = Bk[3*B+j] * Dk[3*A+i]
                        # Term 3
                        term3 = Bk[3*A+i] * Dk[3*B+j]
                        #
                        Ck[3*A+i,3*B+j] = (term1 - term2 - term3) / dom / don
        return np.matrix(Ck)
#-----------------------------------------------#
def wilson_lbend(m,o,n,MON,natoms):
    M,O,N = MON
    #------------------------------------------#
    # Redefine new Cartesian coordinate system #
    #------------------------------------------#
    # Define new z axis
    new_z = n - m
    new_z = new_z / np.linalg.norm(new_z)

    # Define new y axis
    ref  = np.array([+1.0,+1.0,+1.0])
    prod = np.dot(new_z,ref) / np.linalg.norm(new_z) / np.linalg.norm(ref)
    if abs(prod) == 0.0:
       ref = np.array([+1.0,-1.0,0.0])
    new_y = np.cross(ref,new_z)
    new_y = new_y / np.linalg.norm(new_y)

    # Define new x axis
    new_x = np.cross(new_y,new_z)
    new_x = new_x / np.linalg.norm(new_x)

    # Define rotation matrix such as new_r = R * old_r,
    # considering new_r and old_r as column vectors
    R = np.matrix([new_x,new_y,new_z])
    Tbar = np.zeros( (9,9)  )
    Tbar[0:3,0:3] = R
    Tbar[3:6,3:6] = R
    Tbar[6:9,6:9] = R

    # Define position vectors with this new matrix
    m = R * np.matrix(m).transpose(); m = np.array( m.transpose().tolist()[0] )
    o = R * np.matrix(o).transpose(); o = np.array( o.transpose().tolist()[0] )
    n = R * np.matrix(n).transpose(); n = np.array( n.transpose().tolist()[0] )

    #-----------------------------------------------#
    # Get row of B matrix and 2D matrix of C tensor #
    #-----------------------------------------------#
    B_rows , C_matrices = [] , []
    for k in [0,1]:
        Bk, Dk = wilson_auxlinB(m,o,n,k)
        Ck = wilson_auxlinC(m,o,n,k,Bk,Dk)
        # In old Cartesian coord. systems
        Bk = np.matrix(Bk) * Tbar.transpose()
        Ck = Tbar * Ck * Tbar.transpose()
        # Append data
        B_rows.append(np.array(Bk.tolist()[0]))
        C_matrices.append(Ck)

    #------------------------------------------#
    # Consider all the atoms to create B and C #
    #------------------------------------------#
    final_B = []
    for Bk in B_rows:
        fBk = [0.0 for idx in range(3*natoms)]
        for a in [M,O,N]:
            for i in [0,1,2]:
                if a == M: fBk[3*M+i] = Bk[0+i]
                if a == O: fBk[3*O+i] = Bk[3+i]
                if a == N: fBk[3*N+i] = Bk[6+i]
        final_B.append(fBk)

    final_C = []
    for Ck in C_matrices:
        #fCk = np.zeros( (3*natoms,3*natoms) )
        fCk = [ [0.0 for idx1 in range(3*natoms)] for idx2 in range(3*natoms)]
        for a in [M,O,N]:
            for i in [0,1,2]:
                # Select row in C
                if a == M: row = 0+i
                if a == O: row = 3+i
                if a == N: row = 6+i
                for b in [M,O,N]:
                    for j in [0,1,2]:
                        # Select col in C
                        if b == M: col = 0+j
                        if b == O: col = 3+j
                        if b == N: col = 6+j
                        fCk[3*a+i][3*b+j] = Ck[row,col]
        final_C.append(fCk)

    return final_B, final_C
#-----------------------------------------------#
def wilson_torsion(bond_vectors,ijkl,natoms):
    '''
    Returns the row of the B matrix associated to the
    i-j-k-l torsion and the corresponding C matrix.
    Check ref [1] and [5].
    '''

    # Using nomenclature of reference [5]
    m, o, p, n = ijkl
    u, lambda_u = bond_vectors[(o,m)]
    v, lambda_v = bond_vectors[(p,n)]
    w, lambda_w = bond_vectors[(o,p)]

    uxw = np.cross(u,w)
    vxw = np.cross(v,w)

    cosPhi_u = np.dot(u,w) / np.linalg.norm(u) / np.linalg.norm(w)
    sinPhi_u = np.sqrt(1.0 - cosPhi_u**2)
    cosPhi_v = -np.dot(v,w) / np.linalg.norm(v) / np.linalg.norm(w)
    sinPhi_v = np.sqrt(1.0 - cosPhi_v**2)

    # Get internal coordinate: dihedral angle
    cosq = np.dot(uxw,vxw) / sinPhi_u / sinPhi_v
    if   abs(cosq - 1.0) < EPS_SCX:
       cosq = +1.0; q = 0.0
    elif abs(cosq + 1.0) < EPS_SCX:
       cosq = -1.0; q = np.pi
    else: q = np.arccos(cosq)

    # Calculating 1st derivatives: row of B matrix #
    B_row = [0.0 for idx in range(3*natoms)]
    for a in [m,o,p,n]:
        # Get zeta values
        if a == m: zeta_amo = +1.0; zeta_apn =  0.0; zeta_aop =  0.0
        if a == o: zeta_amo = -1.0; zeta_apn =  0.0; zeta_aop = +1.0
        if a == p: zeta_amo =  0.0; zeta_apn = +1.0; zeta_aop = -1.0
        if a == n: zeta_amo =  0.0; zeta_apn = -1.0; zeta_aop =  0.0
        for i in [0,1,2]:
            # Get B element
            dq_dai =  zeta_amo * uxw[i] / lambda_u / sinPhi_u / sinPhi_u + \
                      zeta_apn * vxw[i] / lambda_v / sinPhi_v / sinPhi_v + \
                      zeta_aop * uxw[i] * cosPhi_u / lambda_w / sinPhi_u / sinPhi_u + \
                      zeta_aop * vxw[i] * cosPhi_v / lambda_w / sinPhi_v / sinPhi_v
            B_row[3*a+i] = dq_dai
    # Calculating 2nd derivatives: 2D matrix of C tensor #
    #C_matrix = np.zeros( (3*natoms,3*natoms) )
    C_matrix = [ [0.0 for idx1 in range(3*natoms)] for idx2 in range(3*natoms)]
    for a in [m,o,p,n]:
      for i in [0,1,2]:
        for b in [m,o,p,n]:
          for j in [0,1,2]:
            if C_matrix[3*a+i][3*b+j] != 0.0: continue
            # Define all delta and zeta values
            if a == m: zeta_amo = +1.0; zeta_anp =  0.0; zeta_apo =  0.0; zeta_aop =  0.0; zeta_ano =  0.0
            if a == o: zeta_amo = -1.0; zeta_anp =  0.0; zeta_apo = -1.0; zeta_aop = +1.0; zeta_ano = -1.0
            if a == p: zeta_amo =  0.0; zeta_anp = -1.0; zeta_apo = +1.0; zeta_aop = -1.0; zeta_ano =  0.0
            if a == n: zeta_amo =  0.0; zeta_anp = +1.0; zeta_apo =  0.0; zeta_aop =  0.0; zeta_ano = +1.0

            if b == m: zeta_bom = -1.0; zeta_bmo = +1.0; zeta_bnp =  0.0; zeta_bop =  0.0; zeta_bpo =  0.0
            if b == o: zeta_bom = +1.0; zeta_bmo = -1.0; zeta_bnp =  0.0; zeta_bop = +1.0; zeta_bpo = -1.0
            if b == p: zeta_bom =  0.0; zeta_bmo =  0.0; zeta_bnp = -1.0; zeta_bop = -1.0; zeta_bpo = +1.0
            if b == n: zeta_bom =  0.0; zeta_bmo =  0.0; zeta_bnp = +1.0; zeta_bop =  0.0; zeta_bpo =  0.0
            zeta_bpn = - zeta_bnp

            if a == b: delta_ab = 1.0
            else:      delta_ab = 0.0
            # Get second derivative
            t01 = uxw[i]*(w[j]*cosPhi_u-u[j])/(lambda_u**2)/(sinPhi_u**4)
            t02 = uxw[j]*(w[i]*cosPhi_u-u[i])/(lambda_u**2)/(sinPhi_u**4)
            t03 = vxw[i]*(w[j]*cosPhi_v+v[j])/(lambda_v**2)/(sinPhi_v**4)
            t04 = vxw[j]*(w[i]*cosPhi_v+v[i])/(lambda_v**2)/(sinPhi_v**4)
            t05 = uxw[i]*(w[j]-2*u[j]*cosPhi_u+w[j]*cosPhi_u**2)/(2*lambda_u*lambda_w*sinPhi_u**4)
            t06 = uxw[j]*(w[i]-2*u[i]*cosPhi_u+w[i]*cosPhi_u**2)/(2*lambda_u*lambda_w*sinPhi_u**4)
            t07 = vxw[i]*(w[j]+2*v[j]*cosPhi_v+w[j]*cosPhi_v**2)/(2*lambda_v*lambda_w*sinPhi_v**4)
            t08 = vxw[j]*(w[i]+2*v[i]*cosPhi_v+w[i]*cosPhi_v**2)/(2*lambda_v*lambda_w*sinPhi_v**4)
            t09 = uxw[i]*(u[j]+u[j]*cosPhi_u**2-3*w[j]*cosPhi_u+w[j]*cosPhi_u**3) / (2*lambda_w**2*sinPhi_u**4)
            t10 = uxw[j]*(u[i]+u[i]*cosPhi_u**2-3*w[i]*cosPhi_u+w[i]*cosPhi_u**3) / (2*lambda_w**2*sinPhi_u**4)
            t11 = vxw[i]*(v[j]+v[j]*cosPhi_v**2+3*w[j]*cosPhi_v-w[j]*cosPhi_v**3) / (2*lambda_w**2*sinPhi_v**4)
            t12 = vxw[j]*(v[i]+v[i]*cosPhi_v**2+3*w[i]*cosPhi_v-w[i]*cosPhi_v**3) / (2*lambda_w**2*sinPhi_v**4)
            if i != j and a != b:
               k = [0,1,2]; k.remove(i); k.remove(j); k = k[0]
               t13 = (w[k]*cosPhi_u-u[k]) / lambda_u / lambda_w / sinPhi_u**2
               t14 = (w[k]*cosPhi_v+v[k]) / lambda_v / lambda_w / sinPhi_v**2
            else:
               t13 = 0.0
               t14 = 0.0
            dr_daidbj = zeta_amo*zeta_bmo*(t01 + t02) + \
                        zeta_anp*zeta_bnp*(t03 + t04) + \
                       (zeta_amo*zeta_bop+zeta_apo*zeta_bom)*(t05 + t06) + \
                       (zeta_anp*zeta_bpo+zeta_aop*zeta_bpn)*(t07 + t08) +\
                        zeta_aop*zeta_bpo*(t09 + t10) + \
                        zeta_aop*zeta_bop*(t11 + t12) + \
                       (zeta_amo*zeta_bop+zeta_apo*zeta_bom)*(1-delta_ab)*(j-i) / (-2.)**(abs(j-i))*t13 +\
                       (zeta_anp*zeta_bpo+zeta_aop*zeta_bpn)*(1-delta_ab)*(j-i) / (-2.)**(abs(j-i))*t14
            # Save data in both positions
            C_matrix[3*a+i][3*b+j] = dr_daidbj
            C_matrix[3*b+j][3*a+i] = dr_daidbj
    return [B_row], [C_matrix]
#-----------------------------------------------#
def wilson_getBC(xcc,all_ics):
    numics = count_ics(all_ics)
    ics_st,ics_ab,ics_lb,ics_it,ics_pt = unmerge_ics(all_ics)
    # unit bond vectors
    natoms = fncs.howmanyatoms(xcc)
    bvecs  = wilson_bvecs(xcc)
    # initialize matrices
    wilsonB = [[0.0 for row in range(numics)] for col in range(3*natoms)]
    wilsonB, wilsonC = [], []
    # B and C: stretchings
    for atoms in ics_st:
        row_B, matrix_C = wilson_stretch(bvecs,atoms,natoms)
        wilsonB += row_B
        wilsonC += matrix_C
    # B and C: angular bendings
    for atoms in ics_ab:
        row_B, matrix_C = wilson_abend(bvecs,atoms,natoms)
        wilsonB += row_B
        wilsonC += matrix_C
    # B and C: linear  bendings
    for atoms in ics_lb:
        i,j,k = atoms
        r_m = np.array(fncs.xyz(xcc,i))
        r_o = np.array(fncs.xyz(xcc,j))
        r_n = np.array(fncs.xyz(xcc,k))
        row_B, matrix_C = wilson_lbend(r_m,r_o,r_n,atoms,natoms)
        wilsonB += row_B
        wilsonC += matrix_C
    # B and C: torsions (proper and improper)
    for atoms in ics_it+ics_pt:
        row_B, matrix_C = wilson_torsion(bvecs,atoms,natoms)
        wilsonB += row_B
        wilsonC += matrix_C
    return np.matrix(wilsonB), [np.matrix(ll) for ll in wilsonC]
#-----------------------------------------------#
def wilson_getu(masses):
    '''
       masses: N
       u     : 3Nx3N
    '''
    natoms = len(masses)
    u = np.matrix(np.zeros((3*natoms,3*natoms)))
    for at in range(natoms):
        tt = 1.0/masses[at]
        ii,jj,kk = 3*at+0,3*at+1, 3*at+2
        u[ii,ii] = tt
        u[jj,jj] = tt
        u[kk,kk] = tt
    return u
#-----------------------------------------------#
def wilson_getG(u,B):
    '''
    B     : Fx3N
    G     : FxF
    -------------
    G L  = L Lambda
    G    = L Lambda L^T
    G^-1 = L^T^-1 Lambda^-1 L^-1 = L Lambda^-1 L^T
    -------------
    notice that
       L L^T = L^T L = I
    '''
    # Get G
    G = B*u*B.transpose()
    # Get also inverse
    #Ginv = np.linalg.inv(G) # not valid if singular
    # a) eigenvalues, eigenvectors
    Lambda, L = np.linalg.eigh(G)
    # b) reorder
    idxs = sorted([(abs(li),idx) for idx, li in enumerate(Lambda)])
    idxs = [idx for li,idx in idxs]
    Lambda = [Lambda[idx] for idx in idxs]
    L      = L[:,idxs]
    # Generalized inverse (not inverse of zero eigenvalues)
    G    = L * np.diag(Lambda) * L.transpose()
    Lambda_inv = [1.0/li if li > EPS_GIV else li for li in Lambda]
    Ginv = L * np.diag(Lambda_inv) * L.transpose()
    # return
    return G, Ginv
#-----------------------------------------------#
def wilson_gf_internal(u,B,C,Ginv,gcc,Fcc):
    nIC, n3N = B.shape
    if gcc.shape != (n3N,  1): exit("Wrong shape for gcc")
    if Fcc.shape != (n3N,n3N): exit("Wrong shape for Fcc")
    A = u * B.transpose() * Ginv
    g = A.transpose() * gcc
    f = A.transpose() * Fcc * A
    for i in range(nIC):
        f -=  float(g[i]) * A.transpose() * np.matrix(C[i]) * A
    return g, f, A
#-----------------------------------------------#
def wilson_gf_nonred(G,Ginv,g,f):
    nIC, nIC = G.shape
    if g.shape != (nIC,  1): exit("Wrong shape for g")
    if f.shape != (nIC,nIC): exit("Wrong shape for f")
    P   = G * Ginv
    gnr = P * g
    fnr = P * f * P
    return gnr, fnr
#-----------------------------------------------#
def wilson_prj_rc(gnr,fnr,G,nics):
    ''' project out the reaction coordinate'''
    I   = np.identity(nics)
    p   = gnr * gnr.transpose() / (gnr.transpose() * G * gnr)
    fnr = (I - p*G) * fnr * (I-G*p)
    return fnr
#-----------------------------------------------#
def wilson_evecsincart(L,G,A,masses):
    ''' evectors in cartesian (mass-scaled)'''
    nat = len(masses)
    n3N, nIC = A.shape
    # masses 3N
    m3N = [masses[at] for at in range(nat) for i in range(3)]
    # single value decomposition for L^-1
    u, s, v = np.linalg.svd(L,full_matrices=True,compute_uv=True)
    s_inv = np.diag([1.0/s_i if s_i > EPS_SVD else s_i for s_i in s])
    L_inv = v.transpose() * s_inv.transpose() * u.transpose()
    # Get C matrix
    C = L_inv * G * L_inv.transpose()
    # Get W
    W = np.zeros( C.shape,dtype=complex )
   #for i in range(len(C)): W[i,i] = cmath.sqrt(C[i,i])
    for i in range(len(C)): W[i,i] = np.sqrt(C[i,i])
    # Get chi
    chi = A * L * W
    # Get normal-mode eigenvectors in mass-scaled cartesian
    Lcc  = np.zeros( (3*nat,nIC) )
    for j in range(nIC):
        cocient = sum([m3N[k]*chi[k,j]**2 for k in range(n3N)])
        for i in range(n3N):
            Lcc[i,j] = (np.sqrt(m3N[i])*chi[i,j]/cocient).real
    return Lcc
#-----------------------------------------------#
def calc_icfreqs(Fcc,masses,xcc,gcc,all_ics,bool_prc=False):
    '''
    As described in J. Phys. Chem. A 1998, 102, 242-247
    bool_prc: project reaction coordinate?
    '''
    # some dimensions
    nat = len(masses)
    n3N = 3 * nat
    nIC = count_ics(all_ics)
    # number of vibrational degrees of freedom
    linear = fncs.islinear(xcc)
    if linear: nvdof = n3N - 5
    else     : nvdof = n3N - 6
    if bool_prc: nvdof -= 1
    # the reduced mass (1 amu)
    mu = 1.0/AMU
    # correct shapes:  gcc (3Nx1), Fcc (3Nx3N)
    if bool_prc: gcc = np.matrix(gcc).transpose()
    else       : gcc = np.matrix(np.zeros(n3N)).transpose()
    if len(Fcc) != 3*nat: Fcc = fncs.lowt2matrix(Fcc)
    Fcc = np.matrix(Fcc)
    # 1. Calculate B matrix and C^i tensor
    B, C = wilson_getBC(xcc,all_ics)
    # 2. Calculate G (of h in other paper) and Ginv
    u       = wilson_getu(masses)
    G, Ginv = wilson_getG(u,B)
    # 3. Calculate gradient and Hessian in rics
    g,f,A = wilson_gf_internal(u,B,C,Ginv,gcc,Fcc)
    # 4.1 Project to nrics
    gnr, fnr = wilson_gf_nonred(G,Ginv,g,f)
    # 4.2 Project reaction coordinate
    if bool_prc: fnr = wilson_prj_rc(gnr,fnr,G,nIC)
    # 5.1 Eigenvalues and eigenvectors
    Lambda, L = np.linalg.eig(G*fnr)
    # 5.2 Remove imaginary part and reorder
    Lambda  = [mu*li.real for li in Lambda]
    icfreqs = [fncs.eval2afreq(li,mu) for li in Lambda]
    # 6.1 Transform eigenvecs to mass-scaled Cartesians
    L = wilson_evecsincart(L,G,A,masses)
    # 6.2 save as a list of eigenvectors
    nr,nc = L.shape
    evecs = [L[:,idx].tolist() for idx in range(nc)]
    # 7.1 indices by abs value
    idxs    = sorted([(abs(fq),fq,idx) for idx, fq in enumerate(icfreqs)])
    idxs.reverse()
    # 7.2 keep only the biggest ones (up to the number of degrees of freedom)
    idxs = idxs[:nvdof]
    # 7.3 sort again, imaginary first, then from small to big
    idxs = sorted([(fq,idx) for fqabs,fq,idx in idxs])
    # 7.4 Remove the zero eigenvalues (freqs < 0.1 cm^-1)
    idxs = [(fq,idx) for fq,idx in idxs if abs(fncs.afreq2cm(fq)) > EPS_ICF ]
    # 7.5 Prepare lists
    idxs    = [idx for fq,idx in idxs]
    Lambda  = [Lambda[idx]                   for idx in idxs]
    icfreqs = [icfreqs[idx]                  for idx in idxs]
    evecs   = [evecs[idx]                    for idx in idxs]
    # return data
    return icfreqs, Lambda, evecs
#-----------------------------------------------#
def nonredundant(xcc,masses,gcc,Fcc,all_ics,ccfreqs,unremov=[],ncycles=None,extra=0):
    '''
    unremov: a list with those that cannot be removed
    '''
    mu  = 1.0 / AMU
    nvdof = len(ccfreqs)+extra
    # depure
    count = -1
    original = list(all_ics)
    while True:
          count += 1
          # stop?
          if ncycles is not None and count == ncycles: break
          # select one ic
          target = random.choice(all_ics)
          if target in unremov: continue
          # remove it
          all_ics.remove(target)
          # calculate
          try:
             icfreqs = calc_icfreqs(Fcc,masses,xcc,gcc,all_ics)[0]
             # compare
             same = fncs.same_freqs(ccfreqs,icfreqs)
          except:
             same = False
          # not equal?
          if not same:
              all_ics += [target]
              unremov += [target]
              continue
          # update
          if count_ics(all_ics) <= nvdof : break
    # keep old order
    all_ics = [ic for ic in original if ic in all_ics]
    return all_ics, same
#-----------------------------------------------#
def nonredundant_gtsfiles(gtsfiles,all_ics,unremov=[],ncycles=None,extra=0):
    mu  = 1.0 / AMU
    # depure
    count = -1
    original = list(all_ics)
    while True:
          count += 1
          # stop?
          if ncycles is not None:
             if count == ncycles: break
          # remove one
          target = random.choice(all_ics)
          all_ics.remove(target)
          ic_type, ic_atoms = target
          if ic_atoms in unremov:
              all_ics += [target]
              continue
          # calculate
          for gts in gtsfiles:
              xcc,atonums,ch,mtp,E,gcc,Fcc,masses,pgroup,rotsigma,freq_list = read_gtsfile(gts)
              symbols = fncs.get_symbols(atonums)
              ccfreqs = fncs.calc_ccfreqs(Fcc,masses,xcc)[0]
              icfreqs = calc_icfreqs(Fcc,masses,xcc,gcc,all_ics)[0]
              nvdof   = len(ccfreqs)+extra
             #if count_ics(all_ics) < nvdof:
             #   print "Number of internal coordinates < vibr. degrees of freedom!!"
             #   print "  %i vs %i"%(count_ics(all_ics),nvdof)
             #   exit()
              # compare
              same = fncs.same_freqs(ccfreqs,icfreqs)
              if not same:
                  all_ics += [target]
                  unremov.append(ic_atoms)
                  break
          if not same: continue
          # update
          if count_ics(all_ics) <= nvdof : break
    # keep old order
    all_ics = [ic for ic in original if ic in all_ics]
    return all_ics, same
#-----------------------------------------------#




if __name__ == '__main__': fncs.print_funcs_in_module()
