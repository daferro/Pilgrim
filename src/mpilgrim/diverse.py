#!/usr/bin/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.diverse    |
| Last Update:  2019/01/15 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*

This module contains different
functions used in Pilgrim
'''

#---------------------------------------------------------------#
import os
import names             as     PN
import rdwr              as     RW
import common.fncs       as     fncs
import common.Exceptions as     Exc
from   common.files      import read_gtsfile
from   common.internal   import ic2string
from   common.Molecule   import Molecule
from   common.criteria   import EPS_MEPS
from   common.steepdesc  import sorted_points
#---------------------------------------------------------------#


#===============================================================#
# Reading of all Pilgrim input files at once                    #
#===============================================================#
def get_input_data():
    '''
    * reads all Pilgrim input files
    * it also returns a status list where:
        * -1 : file does not exists
        *  0 : file exists but it is empty
        *  1 : file exists and contains data
    '''
    # read data
    idata  = []
    idata += [ RW.read_ctc()   ] # idx+1 = 1
    idata += [ RW.read_temp()  ] # idx+1 = 2
    idata += [ RW.read_path()  ] # idx+1 = 3
    idata += [ RW.read_tes()   ] # idx+1 = 4
    idata += [ RW.read_chem()  ] # idx+1 = 5
    idata += [ RW.read_kmc()   ] # idx+1 = 6
    idata += [ RW.read_dlevel()] # idx+1 = 7
    # string with status
    string  = ""
    string += "Status of input files\n"
    string += "\n"
    string += "   ---------------------------------\n"
    string += "          input file       | status \n" 
    string += "   ---------------------------------\n"
    idx = 0
    for data,(fname,status) in idata: 
        idx += 1
        string += "    #%i : %-17s |   %2i   \n"%(idx,fname,status)
    string += "   ---------------------------------\n"
    string += "   status = -1 ==> file does not exist\n"
    string += "   status =  0 ==> file exists but it is empty\n"
    string += "   status =  1 ==> file exists and contains data\n"
    # split data
    the_data   = [data   for data,(fname,status) in idata]
    the_status = [status for data,(fname,status) in idata]
    # return data and string
    return the_data, the_status, string
#---------------------------------------------------------------#
def dlevel_to_files(dlevel):
    # files according to case
    dof      = PN.get_dof(dlevel)
    plotfile = PN.get_plf(dlevel)
    hlf      = PN.get_hlf() # high-level file
    # print information
    string  = "Files of importance:\n"
    string += "\n"
    if dlevel     : string += "  - dlevel --> yes\n"
    else          : string += "  - dlevel --> no\n"
    string += "\n"
    string += "  * file for data       storing: '%s'\n"%dof
    string += "  * file for high-level storing: '%s'\n"%hlf
    string += "  * file for plot       storing: '%s'\n"%plotfile
    string += "\n"
    return (dof,hlf,plotfile),string
#===============================================================#


#===============================================================#
# Functions where some things are checked                       #
#===============================================================#
def ffchecking(mustexist=[],tocreate=[]):
    '''
    * checks the existence of folders in 'mustexist'
    * creates folders in 'tocreate'
    '''
    # folders that must exist
    for folder in mustexist:
        if os.path.exists(folder): continue
        return -1
    # folders to create
    for folder in tocreate:
        if not os.path.exists(folder):
           os.mkdir(folder)
    return 0
#---------------------------------------------------------------#
def status_check(lstatus,indices):
    '''
    * checks the status in 'lstatus' associated
      with the indices in filesIDX
    * idx in indices, from 1 to 8
    '''
    final_status = 0
    for idx in indices:
        if lstatus[idx-1] != 1:
           print "     - status of input file #%i is not 1\n"%(idx)
           final_status = -1
    return final_status
#===============================================================#


#===============================================================#
# Functions to do some useful things                            #
#===============================================================#
def pgroup2weight(pgroup):
    #if pgroup.lower() == "c1": return 2
    #if pgroup.lower() == "cs": return 1
    return 1
#---------------------------------------------------------------#
def calc_fwdir(target):
    '''
    * calculates which internal coordinate changes the most
      in the transition state
    * it also gives the sign of the variation in the forward
      direction of the MEP
    '''
    # list of gts files associated with target
    gtsfiles = [PN.DIR1+gts for gts in os.listdir(PN.DIR1)\
                if gts.endswith(".gts") and gts.startswith(target)]
    # Generate Molecule from gts file
    molecule = Molecule()
    molecule.set_from_gts(gtsfiles[0])
    # setup (with frequency calculation)
    molecule.setup()
    ic, fwsign = molecule.get_imag_main_dir()
    # return data
    return (ic2string(ic),fwsign)
#---------------------------------------------------------------#
def get_reaction_energies(TS,dchem,dof):
    '''
    * checks the reactions which involved the target transition
      state (TS) in order to extract V0 and V1 for reactants
      and products 
    '''
    # initialize variables
    reaction = None
    Eref     = None
    V0R      = None
    V0P      = None
    V1R      = None
    V1P      = None
    GibbsR   = None
    # select reaction
    ctc, itc = PN.name2data(TS)
    for rname in dchem.keys():
        Rs, ts, Ps = dchem[rname]
        ctc2, itc2 = PN.name2data(ts)
        if ctc == ctc2:
           reaction = rname
           if itc == itc2: break
    # no reaction?
    if reaction is None: return reaction, V0R, V0P, V1R, V1P, GibbsR
    # read dof
    dall = RW.read_alldata(dof)[0]
    # get energy from reaction
    Rs = dchem[reaction][0]
    Ps = dchem[reaction][2]
    # reactants
    if len(Rs) != 0:
        V0R, V1R = 0.0, 0.0
        for R in Rs:
            ctc, itc = PN.name2data(R)
            if itc is None: key = PN.struckey(ctc,"msho")
            else          : key = R
            data = dall["pfn"].get(key,None)
            if data is None:
               V0R, V1R = None, None
               break
            V0, V1, pfns = data
            V0R += V0
            V1R += V1
            # Gibbs energy
            if key in dall["gibbs"].keys():
               gibbs = dall["gibbs"][key]
               if GibbsR is None: GibbsR = gibbs
               else             : GibbsR = [ gi+gj for gi,gj in zip(GibbsR,gibbs)]
    # products
    if len(Ps) != 0:
        V0P, V1P = 0.0, 0.0
        for P in Ps:
            ctc, itc = PN.name2data(P)
            if itc is None: key = PN.struckey(ctc,"msho")
            else          : key = P
            data = dall["pfn"].get(key,None)
            if data is None:
               V0P, V1P = None, None
               break
            V0, V1, pfns = data
            V0P += V0
            V1P += V1
    # Output
    return reaction, V0R, V0P, V1R, V1P, GibbsR
#---------------------------------------------------------------#
def find_label_in_rst(s_i,drst):
    # initialize output label
    label  = None
    svalue = None
    # get all points
    allpoints = sorted_points(drst,hess=False)
    # as_i is indeed a label!
    if s_i in allpoints: return s_i, drst[s_i][0]
    # s_i is not a float number
    try   : s_i = float(s_i)
    except: return label, s_i
    # go one by one
    for l_j in allpoints:
        s_j  = drst[l_j][0]
        if abs(s_j-s_i) < EPS_MEPS:
           label  = l_j
           svalue = s_j
           break
    return label, svalue
#===============================================================#


