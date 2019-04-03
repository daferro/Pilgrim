#!/usr/bin/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.opt_kmc    |
| Last Update:  2019/01/21 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*

Module for the --kmc option
of Pilgrim
'''

#--------------------------------------------------#
import time
import os
import sys
#--------------------------------------------------#
import numpy             as     np
import rdwr              as     RW
import names             as     PN
import strings           as     PS
#--------------------------------------------------#
from   common.files      import write_file
#--------------------------------------------------#
from   common.fit2anarc  import anarc1
from   common.fit2anarc  import anarc2
from   common.fit2anarc  import anarc3
from   common.fit2anarc  import anarc4
#--------------------------------------------------#
from   common.fncs       import print_string
#--------------------------------------------------#
from   common.kineticsMC import kmc
#--------------------------------------------------#
import common.Exceptions as Exc
#--------------------------------------------------#
from   common.Logger     import Logger
#--------------------------------------------------#
from   common.physcons   import ML
from   common.physcons   import SECOND
from   common.physcons   import NA
#--------------------------------------------------#
from   plotting          import manage_data_for_plot_kmc
from   plotting          import write_plotfile
from   diverse           import get_input_data
from   diverse           import ffchecking
from   diverse           import status_check
#--------------------------------------------------#



#---------------------------------------------------------------#
def get_ratecons(rcs,dchem,dall,idx,temp):
    drcons    = dall.get("rcons",{})
    processes = []
    for key,(rctype,weight,coefs) in sorted(rcs.items()):
        rctype = rctype.lower()
        # reaction name and direction
        if "." in key: rcname, direction = key.split(".")
        else         : rcname, direction = key, "both"
        # elements in reaction
        Rs,TS,Ps = dchem[rcname]
        nR,nP = len(Rs),len(Ps)
        # read/calculate rate constant
        if "analytic" not in rctype:
           key_fw = "%s.%s.%s"%(rctype,rcname,"fw")
           key_bw = "%s.%s.%s"%(rctype,rcname,"bw")
           # get rate constants
           kfw = drcons.get(key_fw,None)
           kbw = drcons.get(key_bw,None)
           if kfw is not None: kfw = kfw[idx]
           if kbw is not None: kbw = kbw[idx]
           # any non-desired rate constants
           if direction == "fw": kbw = None
           if direction == "bw": kfw = None
           # none
           if kbw is None and kfw is None:
              exception =Exc.NoRateCons(Exception)
              exception._var = (rctype,key)
              raise exception
        else:
           if   rctype.lower() == "analytic1": k = anarc1(temp,*coefs)
           elif rctype.lower() == "analytic2": k = anarc2(temp,*coefs)
           elif rctype.lower() == "analytic3": k = anarc3(temp,*coefs)
           elif rctype.lower() == "analytic4": k = anarc4(temp,*coefs)
           else                              : k = None
           # save data properly
           if   direction in ["fw","both"]: kfw, kbw = k   , None
           elif direction in ["bw"       ]: kfw, kbw = None, k
           else                           : kfw, kbw = None, None
           # in atomic units
           hunitsFW = ML**(nR-1.0) / SECOND
           hunitsBW = ML**(nP-1.0) / SECOND
           if kfw is not None: kfw /= hunitsFW
           if kbw is not None: kbw /= hunitsBW
        # ignore reactions giving rise to bimolecular products
        if len(Ps) > 1 and direction != "bw": kbw = None
        # save in processes
        if kfw is not None: processes.append( (Rs,Ps,weight*kfw) )
        if kbw is not None: processes.append( (Ps,Rs,weight*kbw) )
    return processes
#---------------------------------------------------------------#
def kmc_savedata(kmcof,data):
    # Molecules and same length
    molecules = sorted(data[0][2].keys())
    len1 = max([len(string) for string in molecules+["time","-"*10]])
    ss   = "%%%is"%len1
    # string for gnuplot data file
    string = ""
    string += "# Time in ps; populations in molecules\n"
    string += "\n\n"
    # string for gnuplot string
    gpstr  = "datafile='%s'\n"%kmcof
    gpstr += "plot\\\n"
    # data for each temperature
    for ii,(stemp,xvalues,yvalues) in enumerate(data):
         # add line in gnuplot script
         gpstr += "    datafile  index %i  using 1:2 with points,\\\n"%ii
         # add data in gnuplot file
         string += "# T = %s K\n"%stemp
         string += "#"+ss%"time" + " ".join([ss%mol for mol in molecules]) + "\n"
         for jj,time in enumerate(xvalues):
             time = 1e12*SECOND*time
             string += "  %10.3E "%time
             for molecule in molecules:
                 pop = yvalues[molecule][jj]/NA
                 string += " %10.3E "%pop
             string += "\n"
         string += "\n\n"
    # write files
    write_file(kmcof,string)
    write_file("plotkmc.gp",gpstr)
#---------------------------------------------------------------#



def main(idata,status,case):

    stat2check = [2,5,6]
    mustexist  = []
    tocreate   = [PN.DIR3,PN.DIR6]
    #-------------------------------------------------------#
    # Read Pilgrim input files, check file/folder status    #
    # and expand tuple 'case'                               #
    #-------------------------------------------------------#
    # expand data
    (dctc,dimasses), ltemp, dpath, (dtesLL,dtesHL), dchem, tkmc, ddlevel = idata
    # status ok?
    fstatus = status_check(status,stat2check)
    if fstatus == -1: exit()
    # existency of folders
    fstatus = ffchecking(mustexist,tocreate)
    if fstatus == -1: exit()
    # expand case
    (dof,hlf,plotfile),dlevel,software = case
    #-------------------------------------------------------#


    #---------------------------------#
    # files with data and output file #
    #---------------------------------#
    pof = PN.get_pof(dlevel,"kmc")
    print "   Pilgrim output file: %s"%pof
    print

    sys.stdout = Logger(pof,"w",True)

    # expand KMC tuple
    ipops,rcs,psteps,volume,timeunits,excess = tkmc
    valid_tu = "fs,ps,mcs,ms,s,min,hr".split()
    if timeunits not in valid_tu: valid_tu = "ps"

    # reactivo limitante
    try   : POP0 = min([ipop for species,ipop in sorted(ipops.items()) if ipop != 0.0])
    except: POP0 = 0

    print_string(PS.skmc_init(ipops,POP0,excess,rcs,psteps,volume,timeunits,dof),5)

    # continue?
    if POP0 == 0:
        print "   All initial populations are zero..."
        print
        return
    if rcs == {}:
        print "   No reactions considered in %s"%PN.IFILE6
        print
        return

    # read dof
    dall = RW.read_alldata(dof,ltemp)[0]

    # perform KMC for each temperature
    data = []
    fratios = {}
    ftimes  = []
    for idx,temp in enumerate(ltemp):
        # title
        title    = " T = %.3f Kelvin "%temp
        division ="-"*len(title)
        string  = "     "+division+"\n"
        string += "     "+title   +"\n"
        string += "     "+division+"\n"
        print string
        # get rate constants
        processes = get_ratecons(rcs,dchem,dall,idx,temp)
        # print initial information before simulation
        print_string(PS.skmc_processes(temp,processes),9)
        # perform kmc simulation
        xvalues, yvalues = kmc(ipops,processes,excess,volume,psteps)
        fratios["%7.2f"%temp] = {}
        for species in yvalues.keys():
            fratios["%7.2f"%temp][species] = yvalues[species][-1] / POP0
        # print data from simulation
        print_string(PS.skmc_results(xvalues,yvalues,timeunits),9)
        # save data needed for txt and pdf files
        data += [ ("%.2f"%temp,xvalues,yvalues) ]
        ftimes.append(xvalues[-1])

    # print final ratios
    species = sorted(yvalues.keys())
    print_string(PS.skmc_finaltimes(ltemp,ftimes,timeunits),5)
    print_string(PS.skmc_finalratio(ltemp,species,fratios),5)

    # save data for plotting
    if plotfile is not None:
       plotdata = {}
       plotdata.update(manage_data_for_plot_kmc(data,fratios,volume,timeunits))
       write_plotfile(plotfile,plotdata)
    #kmc_savedata(kmcof,data)

#---------------------------------------------------------------#


