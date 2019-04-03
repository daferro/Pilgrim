#!/usr/bin/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.itf_asurf  |
| Last Update:  2019/01/15 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*

InTerFace betweem Pilgrim and an Analytic Surface
'''

#---------------------------------------------------------------#
import os
import common.fncs    as fncs
import common.asurf   as ITF
#---------------------------------------------------------------#


#===============================================================#
def pilgrim_template(ch=0,mtp=1,nproc=1,mem=1):
    '''
    Input for analytic surface
    This calculation IS a single point calculation!
    '''
    return None
#---------------------------------------------------------------#
def pilgrim_spc(xcc,symbols,bhessian,mname,eargs):
    atonums = fncs.symbols2atonums(symbols)
    ch, mtp = None, None
    masses  = None
    # expand extra-args
    if   len(eargs) == 1:
        asurf = eargs[0]
        folder, clean = None, False
    elif len(eargs) == 2:
        asurf, folder = eargs
        clean = False
    elif len(eargs) == 3:
        asurf, folder, clean = eargs
    # create folder
    if folder is not None:
        if not os.path.exists(folder): os.mkdir(folder)
    # names of files
    wname, ifile, ofile, err = ITF.iofiles(mname,folder)
    # write ifile
    ITF.write_input(ifile,xcc)
    # calculate gradient/hessian
    if bhessian: E, gcc, Fcc = ITF.calculate_hessian(ifile,ofile,err,asurf)
    else       : E, gcc, Fcc = ITF.calculate_gcc(ifile,ofile,err,asurf)+[None]
    # Remove files
    if clean:
       files = os.listdir(folder)
       for fff in files: os.remove(folder+fff)
    return xcc, atonums, ch, mtp,  E, gcc, Fcc, masses
#===============================================================#


