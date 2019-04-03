#!/usr/bin/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.names      |
| Last Update:  2019/02/05 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*

This module contains the names of
files and folders used in Pilgrim.

There are also functions to get
such a names.
'''


#===============================================================#
UFOLDER = "UDATA/"                                              #
ANHDIR  = "ANHAR/"                                              #
TMP     = "TMP/"                                                #
#===============================================================#
TMPi   = TMP + "MEPcalcs_%s/"                                   #
TMPHLi = TMP + "HLcalcs_%s/"                                    #
DIR1   = "1-GTS/"                                               #
DIR2   = "2-PLG_DATA/"                                          #
DIR3   = "3-PLG_OUTPUT/"                                        #
DIR4   = "4-PLG_RST/"                                           #
DIR5   = "5-MOLDEN/"                                            #
DIR6   = "6-PLOTFILES/"                                         #
#===============================================================#
IFILE0 = "tracking"                                             #
IFILE1 = "pif.struc"                                            #
IFILE2 = "pif.temp"                                             #
IFILE3 = "pif.path"                                             #
IFILE4 = "pif.calcs"                                            #
IFILE5 = "pif.chem"                                             #
IFILE6 = "pif.kmc"                                              #
IFILE7 = "pif.dlevel"                                           #
#===============================================================#


#===============================================================#
def name2data(name):                                  
    if "." in name: ctc, itc = name.split(".")        
    else          : ctc, itc = name, None             
    return ctc, itc                                   
#---------------------------------------------------------------#
def struckey(ctc,itc=None):
    name = ctc                                        
    if itc is not None: name += "."+itc               
    return name                                       
#===============================================================#


#===============================================================#
def get_gts(ctc,itc):
    return DIR1 + "%s.gts"%struckey(ctc,itc)
#---------------------------------------------------------------#
def get_rst(ctc,itc):
    name = struckey(ctc,itc)
    return DIR4 + "%s.rst"%name
#---------------------------------------------------------------#
def get_dof(dlevel):
    if dlevel: fname = "data.dlevel"
    else     : fname = "data.slevel"
    return DIR2 + fname
#---------------------------------------------------------------#
def get_hlf():
    fname = "highlevel.txt"
    return DIR2 + fname
#---------------------------------------------------------------#
def get_pof(dlevel,case,target=None):
    fname = case
    if target is not None: fname += ".%s"%target
    if dlevel            : fname += ".dlevel.txt"
    else                 : fname += ".slevel.txt"
    return DIR3 + fname
#---------------------------------------------------------------#
def get_gtsmolden(ctc,itc):
    return DIR5 + "sp.%s.molden"%struckey(ctc,itc)
#---------------------------------------------------------------#
def get_rstxyz(ctc,itc):
    name = struckey(ctc,itc)
    return DIR5 + "path.%s.xyz"%name
#---------------------------------------------------------------#
def get_plf(dlevel=False):
    if dlevel: fname = "plots.dlevel.txt"
    else     : fname = "plots.slevel.txt"
    return DIR6 + fname
#===============================================================#


