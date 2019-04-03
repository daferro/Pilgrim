#!/usr/bin/python2.7
'''
*----------------------------------*
| Module name:  asurf              |
| Last Update:  2018/12/05 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
Interface for the Electronic Structure calculation
with an ANALYTIC SURFACE

# How to use analytic surfaces:
surface    : asurf.x
execution  : asurf.x  ifile  ofile
input_info : atomic coordinates [in bohr]
output_info: energy and gradient

# Input example
   0.00000  0.00000  0.00000
   0.30000  0.20000  0.00000
   0.30000 -0.20000  0.00000
   0.10000 -0.10000  0.20000

# Output example
   -134.6756374
   0.50000 -0.00030  0.13400
   0.00000  1.20000  1.24000
   0.00000 -3.10000  0.00000
  -0.08000 -1.10000 -0.10000

'''

import os
import math
import fncs       as fncs
import files      as ff
import Exceptions as Exc
from   criteria import EPS_HESSDX

WAS_SET     = False
WAS_CHECKED = False

#======== EXE files ========#
dEXE = { \
"clhbrpot":"/home/david/Dropbox/PyFerro/tests/clhbrpot/data/clhbrpot.x" \
}
#===========================#

#==========================================================#
def set_EXE(asurf=None):
    global EXE, WAS_SET
    WAS_SET = True
    # Defined in this file
    if (asurf is None) or ('EXE' in globals()): return
    # Get path for analytic surface
    asurf = asurf.strip()
    if asurf in dEXE.keys(): EXE = dEXE[asurf]
    else                   : raise Exc.ExeNotDef(Exception)
#----------------------------------------------------------#
def check_EXE():
    global WAS_CHECKED
    WAS_CHECKED = True
    if ('EXE' not in globals()): raise Exc.ExeNotDef(Exception)
    if not os.path.exists(EXE) : raise Exc.ExeNotFound(Exception)
#----------------------------------------------------------#
def iofiles(name,folder=None):
    '''
    For a given name, it returns the name of all files
    for the calculation
    '''

    if   folder is None          : folder = ""
    elif not folder.endswith("/"): folder = folder+"/"
    else                         : folder = folder
    wname  = folder + name # whole name
    ifile  = wname + ".iasurf"
    ofile  = wname + ".oasurf"
    err    = wname + ".err"
    return wname, ifile, ofile, err
#----------------------------------------------------------#
def read_input(ifile):
    lines = ff.read_file(ifile)
    xcc   = [line.split() for line in lines]
    xcc   = fncs.flatten_llist(xcc)
    xcc   = [float(x_i) for x_i in xcc]
    return xcc
#----------------------------------------------------------#
def write_input(ifile,xcc,f="%+13.9f"):
    string = ""
    natoms = fncs.howmanyatoms(xcc)
    for atom in range(natoms):
        x = f%(fncs.x(xcc,atom))
        y = f%(fncs.y(xcc,atom))
        z = f%(fncs.z(xcc,atom))
        string += "  %s %s %s \n"%(x,y,z)
    ff.write_file(ifile,string)
#----------------------------------------------------------#
def read_output(ofile):
    lines = ff.read_file(ofile)
    E     = float(lines[0])
    gcc   = [line.split() for line in lines[1:]]
    gcc   = fncs.flatten_llist(gcc)
    gcc   = [float(g_i) for g_i in gcc]
    return E, gcc
#----------------------------------------------------------#
def execute(ifile,ofile,err,asurf=None):
    # Execution file?
    if not WAS_SET    : set_EXE(asurf)
    if not WAS_CHECKED: check_EXE()
    # Exception
    exception = Exc.CalcFails(Exception)
    exception._var = ifile
    # Execute analytic surface
    command = "%s %s %s 2>%s"%(EXE,ifile,ofile,err)
    try   : status  = os.system(command)
    except: raise exception
    return status
#----------------------------------------------------------#
def calculate_gcc(ifile,ofile,err,asurf=None):
    '''
    single point calculation (just gradient)
    '''
    # Execution
    execute(ifile,ofile,err,asurf)
    # Read output
    E, gcc = read_output(ofile)
    # Return data
    return [E, gcc]
#----------------------------------------------------------#
def calculate_hessian(ifile,ofile,err,asurf=None):
    '''
    single point calculation (with numerical hessian)
    '''
    # name from ifile
    name = ".".join(ifile.split(".")[:-1])
    # ifile info
    xcc    = read_input(ifile)
    natoms = fncs.howmanyatoms(xcc)
    # Energy and gradient
    E, gcc = calculate_gcc(ifile,ofile,err,asurf)
    # Hessian
    hessian = []
    for coord in range(len(xcc)):
        for idx in range(6):
            #----------------#
            # Left  geometry #
            #----------------#
            xcc_L        = list(xcc)
            xcc_L[coord] = xcc_L[coord] - EPS_HESSDX
            # Left  file
            nameL, ifileL, ofileL, efileL = iofiles(name+".L",folder=None)
            if idx == 0: write_input(ifileL,xcc_L,f="%+13.9f")
            if idx == 1: write_input(ifileL,xcc_L,f="%+13.8f")
            if idx == 2: write_input(ifileL,xcc_L,f="%+13.7f")
            if idx == 3: write_input(ifileL,xcc_L,f="%+13.6f")
            if idx == 4: write_input(ifileL,xcc_L,f="%+13.5f")
            if idx == 5: write_input(ifileL,xcc_L,f="%+13.4f")
            # Left  calculation
            E_L, gcc_L = calculate_gcc(ifileL,ofileL,efileL,asurf)
            #----------------#
            # Right geometry #
            #----------------#
            xcc_R        = list(xcc)
            xcc_R[coord] = xcc_R[coord] + EPS_HESSDX
            # Right file
            nameR, ifileR, ofileR, efileR = iofiles(name+".R",folder=None)
            if idx == 0: write_input(ifileR,xcc_R,f="%+13.9f")
            if idx == 1: write_input(ifileR,xcc_R,f="%+13.8f")
            if idx == 2: write_input(ifileR,xcc_R,f="%+13.7f")
            if idx == 3: write_input(ifileR,xcc_R,f="%+13.6f")
            if idx == 4: write_input(ifileR,xcc_R,f="%+13.5f")
            if idx == 5: write_input(ifileR,xcc_R,f="%+13.4f")
            # Right calculation
            E_R, gcc_R = calculate_gcc(ifileR,ofileR,efileR,asurf)

            #----------------#
            # Data 2 hessian #
            #----------------#
            row = [(gR-gL)/(2.0*EPS_HESSDX) for gR,gL in zip(gcc_R,gcc_L)]

            # Recalculate (any NaN in row?)
            newtry = False
            for rr in row:
                if math.isnan(float(rr)): newtry = True
            if newtry and idx<6: continue

            # Append data
            hessian.append(row)
            break

    # Average hessian
    nrows = len(hessian)
    ncols = len(hessian[0])
    for i in range(1,nrows,1):
        for j in range(i):
           Fij = hessian[i][j]
           Fji = hessian[j][i]
           average = 0.5*(Fij + Fji)
           hessian[i][j] = average
           hessian[j][i] = average
    # Return
    return E, gcc, fncs.matrix2lowt(hessian)
#==========================================================#

