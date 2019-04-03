#!/usr/bn/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.itf        |
| Last Update:  2019/01/15 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*

This module is in charge of returning
the Python function used for the
calculation with selected software.
'''


#---------------------------------------------------------------#
import os
from   common.files      import read_file
from   common.fncs       import clean_lines
from   common.Exceptions import UnknownSoft
#---------------------------------------------------------------#


#===============================================================#
def get_dsoft():
    '''
    This functions is in charge of reading the file
    mesc.txt where the available softwares are listed
    '''
    # defined here:
    if False:
        dsoft = {'clhbrpot': 'mpilgrim.itf_asurf', \
                 'gaussian': 'mpilgrim.itf_gau'  ,\
                 'orca'    : 'mpilgrim.itf_orca'   \
                }
    # read from file
    if True:
       TXTFILE = os.path.dirname(os.path.realpath(__file__))+"/mesc.txt"
       # read file
       lines = read_file(TXTFILE)
       lines = clean_lines(lines,"#",True)
       # get info from lines
       dsoft   = {}
       for line in lines:
           if line == "\n": continue
           software, module = line.split()
           dsoft[software] = module
    # return dictionary
    return dsoft
#---------------------------------------------------------------#
def get_templates(ch,mtp,case="LL"):
    templates = {}
    # get info from file
    dsoft = get_dsoft()
    for software in dsoft.keys():
        module = dsoft[software]
        lib = __import__(module,fromlist=[''])
        globals()["mes"] = lib
        if   case == "LL": string = mes.pilgrim_template(ch,mtp)
        elif case == "HL": string = mes.pilgrim_templateHL(ch,mtp)
        else             : string = ""
        if string is None: string = software+"\n" # for asurf
        templates[software] = string
    return templates
#---------------------------------------------------------------#
def get_spc_fnc(software):
    '''
    Returns function to performs a single-point calculation
    using the corresponding interface
    '''
    # get info from file
    dsoft = get_dsoft()
    # valid software?
    if software not in dsoft.keys():
        raise UnknownSoft(Exception)
    # import module
    module = dsoft[software]
    lib = __import__(module,fromlist=[''])
    globals()["mesc"] = lib
    # call function from module and return data
    return mesc.pilgrim_spc
#===============================================================#


