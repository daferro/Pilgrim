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
| Module name:  orca               |
| Last Update:  2019/04/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
Interface for the Electronic Structure calculation
using ORCA
'''

#==========================================================#
#                  TO BE MODIFIED BY USER                  #
#==========================================================#
#EXE = "/home/david/Software/orca_4_0_1_2/orca"             #
#==========================================================#


#==========================================================#
import os                                                  #
import Exceptions as Exc
from   fncs      import xyz
from   fncs      import clean_lines
from   physcons  import ANGSTROM
from   files     import read_file
from   files     import write_file
#==========================================================#


#==========================================================#
def set_EXE():
    global EXE
    txt = os.path.dirname(os.path.realpath(__file__))+"/paths.txt"
    # Defined in this file
    if 'EXE' in globals():
        return
    # Try to export it from bashrc
    elif "OrcaExe" in os.environ:
        # in .bashrc: export OrcaExe="$MYHOME/Software/orca_4_0_1_2/orca"
        EXE = os.environ["OrcaExe"]
        return
    # Export it from file
    elif os.path.exists(txt):
        lines = read_file(txt)
        lines = clean_lines(lines,"#",True)
        for line in lines:
            if line == "\n": continue
            name, path = line.split()
            if name == "orca":
                EXE = path
                return
    # Not found
    else: raise Exc.ExeNotDef(Exception)
#----------------------------------------------------------#
def check_EXE():
    if not os.path.exists(EXE): raise Exc.ExeNotFound(Exception)
#----------------------------------------------------------#
def execute(ifile,ofile,err,folder=None):
    # Add folder to names
    if folder is not None:
       if not folder.endswith("/"): folder = folder + "/"
       ifile = folder + ifile
       ofile = folder + ofile
       err   = folder + err
    # Execution file?
    set_EXE()
    check_EXE()
    # Exception
    exception = Exc.CalcFails(Exception)
    exception._var = ofile
    # Execute command
    command = "%s %s 1>%s 2>%s"%(EXE,ifile,ofile,err)
    try   : status  = os.system(command)
    except: raise exception
    return status
#----------------------------------------------------------#
def normal_termination(ofile):
    olines = read_file(ofile)
    if len(olines) == 0: return False
    if "ORCA TERMINATED NORMALLY" not in olines[-2]:
        return False
    return True
#----------------------------------------------------------#
def iofiles(name,folder=None):
    '''
    For a given name, it returns the name of all files
    for the Gaussian calculation
    '''

    if   folder is None          : folder = ""
    elif not folder.endswith("/"): folder = folder+"/"
    else                         : folder = folder
    wname  = folder + name # whole name
    ifile  = wname + ".inp"
    ofile  = wname + ".out"
    engrad = wname + ".engrad"
    hess   = wname + ".hess"
    gbw    = wname + ".gbw"
    err    = wname + ".err"
    return wname, ifile, ofile, engrad, hess, gbw, err
#==========================================================#


