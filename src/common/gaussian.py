#!/usr/bin/python2.7
'''
*----------------------------------*
| Module name:  gaussian           |
| Last Update:  2018/12/05 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
Interface for the Electronic Structure calculation
using GAUSSIAN
'''

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
#                  TO BE MODIFIED BY USER                  #
#==========================================================#
#EXE  = "/home/programs/G09_64david/g09/g09"               #
#FCHK = "/home/programs/G09_64david/g09/formchk"           #
#----------------------------------------------------------#
# in .bashrc:                                              #
#  export GauExe="/home/programs/G09_64david/g09/g09"      #
#  export GauFchk="/home/programs/G09_64david/g09/formchk" #
#==========================================================#


#==========================================================#
def set_EXE():
    global EXE
    # Defined in this file
    if 'EXE' in globals():
        return
    # Try to export it from bashrc
    elif "GauExe" in os.environ:
        EXE = os.environ["GauExe"]
        return
    # Not found
    else: raise Exc.ExeNotDef(Exception)
#----------------------------------------------------------#
def set_FCHK():
    global FCHK
    # Defined in this file
    if 'FCHK' in globals():
        return
    # Try to export it from bashrc
    elif "GauFchk" in os.environ:
        FCHK = os.environ["GauFchk"]
        return
    # Not found
    else: raise Exc.ExeNotDef(Exception)
#----------------------------------------------------------#
def check_EXE():
    if not os.path.exists(EXE): raise Exc.ExeNotFound(Exception)
#----------------------------------------------------------#
def check_FCHK():
    if not os.path.exists(FCHK): raise Exc.ExeNotFound(Exception)
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
    # Execute Gaussian
    command = "%s <%s 1>%s 2>%s"%(EXE,ifile,ofile,err)
    try   : status  = os.system(command)
    except: raise exception
    return status
#----------------------------------------------------------#
def normal_termination(ofile):
    lines = read_file(ofile)
    if len(lines) == 0: return False
    lastline = lines[-1].lower()
    if "normal termination" in lastline: return True
    else                               : return False
#----------------------------------------------------------#
def genfchk(chk,fchk,err,folder=None):
    # Add folder to names
    if folder is not None:
       if not folder.endswith("/"): folder = folder + "/"
       chk   = folder + chk
       fchk  = folder + fchk
       err   = folder + err
    # fchk tool?
    set_FCHK()
    check_FCHK()
    # Exception
    exception = Exc.CalcFails(Exception)
    exception._var = chk
    # Execute fchk tool
    command = "%s %s %s 1>%s 2>&1"%(FCHK,chk,fchk,err)
    try   : status = os.system(command)
    except: raise exception
    return status
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
    ifile  = wname + ".gjf"
    ofile  = wname + ".out"
    chk    = wname + ".chk"
    fchk   = wname + ".fchk"
    err    = wname + ".err"
    return wname, ifile, ofile, chk, fchk, err
#==========================================================#


