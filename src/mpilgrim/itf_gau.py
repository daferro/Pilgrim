#!/usr/bin/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.itf_gau    |
| Last Update:  2019/01/15 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*

InTerFace betweem Pilgrim and Gaussian
'''

#---------------------------------------------------------------#
import os                            
import time
import common.Exceptions   as Exc    
import common.fncs         as fncs   
import common.physcons     as pc     
import common.files        as ff     
import common.gaussian     as ITF    
#---------------------------------------------------------------#


#===============================================================#
key1 = "[Pilgrim_geometry]"
key2 = "[Pilgrim_name]" 
key3 = "[Pilgrim_gradhess]"
#---------------------------------------------------------------#
def pilgrim_template(ch=0,mtp=1,nproc=1,mem=1):
    '''
    Input for Gaussian software
    This calculation IS a single point calculation!
    '''
    level   = "hf/sto-3g"
   #level   = "mpwb95/6-31+G(d,p) IOp(3/76=0560004400)"

    string  = "%%nproc=%i   \n"%nproc
    string += "%%mem=%iGB   \n"%mem
    string += "%%chk=%s.chk \n"%key2
    string += "#p %s        \n"%level
    string += "scf=verytight\n"
   #string += "int=ultrafine\n"
    string += "NoSymm       \n"
    string += "%s           \n"%key3
    string += "\n"
    string += "Input file for MEP calculation\n"
    string += "\n"
    string += "%i %i        \n"%(ch,mtp)
    string += "%s           \n"%key1
    string += "\n"
    return string
#---------------------------------------------------------------#
def pilgrim_templateHL(ch=0,mtp=1,nproc=1,mem=1):
    '''
    Input for Gaussian software
    This calculation IS a single point calculation!
    '''
    level   = "hf/6-31G"

    string  = "%%nproc=%i   \n"%nproc
    string += "%%mem=%iGB   \n"%mem
    string += "%%chk=%s.chk \n"%key2
    string += "#p %s        \n"%level
    string += "scf=verytight\n"
    string += "NoSymm       \n"
    string += "\n"
    string += "Input file for MEP calculation\n"
    string += "\n"
    string += "%i %i        \n"%(ch,mtp)
    string += "%s           \n"%key1
    string += "\n"
    return string
#---------------------------------------------------------------#
def pilgrim_spc(xcc,symbols,bhessian,mname,eargs):
    # expand extra-args
    if   len(eargs) == 1:
        spc_template = eargs[0]
        folder, clean = None, False
    elif len(eargs) == 2:
        spc_template, folder = eargs
        clean = False
    elif len(eargs) == 3:
        spc_template, folder, clean = eargs
    # no template??
    if spc_template is None: raise Exc.NoTemplateGiven(Exception)
    # create folder
    if folder is not None:
        if not os.path.exists(folder): os.mkdir(folder)
    # Calculate gradient&hessian or just gradient
    if bhessian: inkey3 = "freq=noraman "
    else       : inkey3 = "force        "
    # names of files
    wname, ifile, ofile, chk, fchk, err = ITF.iofiles(mname,folder)
    # Input
    string_ifile = ""
    for line in spc_template.split("\n"):
        if key1 in line:
           line = ""
           for idx,symbol in enumerate(symbols):
               x,y,z  = fncs.xyz(xcc,idx)
               linedata = (symbol,x*pc.ANGSTROM,y*pc.ANGSTROM,z*pc.ANGSTROM)
               line += "%2s   %+11.6f  %+11.6f  %+11.6f\n"%linedata
        if key2 in line:
           pos  = line.find(key2)
           line = line[0:pos] + wname + line[pos+len(key2):]
        if key3 in line:
           pos  = line.find(key3)
           line = line[0:pos] + inkey3 + line[pos+len(key3):]
        # Add \n to line
        if not line.endswith("\n"): line += "\n"
        string_ifile += line
    ff.write_file(ifile,string_ifile)
    # Send calculation
    status = ITF.execute(ifile,ofile,err)
 #  # Sleep 0.5 seconds, so Gaussian can write the files
 #  time.sleep(0.5)
    # Check output
    exception = Exc.CalcFails(Exception)
    exception._var = ofile
    if not ITF.normal_termination(ofile): raise exception
    # Generate fchk
    status = ITF.genfchk(chk,fchk,err)
    # Read data
    xcc, atonums, ch, mtp, E, gcc, Fcc, masses, clevel = ff.read_fchk(fchk)
    # Remove files
    if clean:
       files = os.listdir(folder)
       files = [fff for fff in files if     fff.startswith(name)   ]
       files = [fff for fff in files if not fff.endswith(".gjf")   ]
       files = [fff for fff in files if not fff.endswith(".out")   ]
       files = [fff for fff in files if not fff.endswith(".chk")   ]
       files = [fff for fff in files if not fff.endswith(".fchk")  ]
       for fff in files: os.remove(folder+fff)
    return xcc, atonums, ch, mtp,  E, gcc, Fcc, masses
#===============================================================#



