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

*-----------------------------------*
| Module name:  mpilgrim.itf_asurf  |
| Last Update:  2019/04/03 (Y/M/D)  |
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


