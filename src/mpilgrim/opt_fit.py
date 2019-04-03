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
| Module name:  mpilgrim.opt_fit    |
| Last Update:  2019/04/03 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*

Module for the --fit option
of Pilgrim
'''

#--------------------------------------------------#
import names   as PN
import rdwr    as RW
import strings as PS
#--------------------------------------------------#
import common.Exceptions as Exc
import common.partfns    as partfns
#--------------------------------------------------#
from   common.fit2anarc  import fit2anarc
from   common.fncs       import print_string
import common.physcons as pc
#--------------------------------------------------#
from   diverse           import ffchecking
from   diverse           import status_check
from   plotting          import manage_data_for_plot_rcons
from   plotting          import write_plotfile
#--------------------------------------------------#


#===============================================================#
def fit2analytic(ltemp,rcons,nR):
    # fit to analytic
    if len(ltemp) > 1:
       human_units = pc.ML**(nR-1.0) / pc.SECOND
       k    = [val*human_units for val in rcons]
       dfit = fit2anarc(ltemp,k,log=True)
       # save data
       return (k,dfit)
    # return data
    return (None,None)
#===============================================================#



#===============================================================#
def main(idata,status,case,targets="*"):

    stat2check = [1,2,5]
    mustexist  = [PN.DIR1]
    tocreate   = [PN.DIR2,PN.DIR3,PN.DIR6]
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


    # no specific target selected
    if "*" in targets or len(targets) == 0: targets = dchem.keys()

    # clean targets
    targets = sorted([target for target in targets if target in dchem.keys()])

    # read dof
    dall = RW.read_alldata(dof,ltemp)[0]
    
    #--------------------------------#
    # Calculations for each reaction #
    #--------------------------------#
    print  "   Analytic expressions:"
    print
    print  "       (1) k = A exp(-B/T)              "
    print  "       (2) k = A*T^n*exp(-B/T)          "
    print  "       (3) k = A*(T/Tr)^n*exp(-B/T)     "
    print  "       (4) k = A*(T/Tr)^n*exp(-B*(T+T0)/(T^2+T0^2))"
    print
    plotdata = {}
    LINES = {}
    for rcname in targets:
        Rs, TS, Ps = dchem[rcname]
        nR = len(Rs)
        nP = len(Ps)
        title = " Reaction Name: %s "%rcname
        print "   "+"-"*len(title)
        print "   "+        title
        print "   "+"-"*len(title)
        print
        for direc,num in [("fw",nR),("bw",nP)]:
            keys  = [key for key in dall["rcons"].keys() if key.split(".")[1:3] == [rcname,direc] ]
            if keys == []: continue
            # print table of forward rate constants
            if direc == "fw":
                print "    * Reaction in FORWARD dir. is '%s --> %s'"%(" + ".join(Rs)," + ".join(Ps))
                print
            # print table of backward rate constants
            if direc == "bw":
                print "    * Reaction in BACKWARD dir. is '%s --> %s'"%(" + ".join(Ps)," + ".join(Rs))
                print
            # harmonic and anharmonic keys
            ho_keys  = [key for key in keys if not key.startswith("anh")]
            anh_keys = [key for key in keys if     key.startswith("anh")]
            for keys,case in [(ho_keys,0),(anh_keys,1)]:
                if keys == []: continue
                rcons = {}
                for key in keys:
                    rctype = key.split(".")[0]
                    RC  = dall["rcons"][key]
                    rcons[rctype] = RC
                d4plot = {}
                print_string(PS.sfit_rcons(ltemp,rcons,num,case),6)
                # string for units
                if num != 1: units = "s^{-1} cc^{%i} molecule^{-%i}"%(num-1,num-1)
                else       : units = "s^{-1}"
                # Fitting process
                print "      Fitting parameters:"
                print
                for rctype in "tst,cvt,tstzct,cvtzct,tstsct,cvtsct".split(","):
                    if case == 0: key = "%s.%s.%s"%(rctype,rcname,direc)
                    else        : key = "anh%s.%s.%s"%(rctype,rcname,direc)
                    if key not in dall["rcons"].keys(): continue
                    rcons    = dall["rcons"][key]
                    # the fitting
                    khu, dfit = fit2analytic(ltemp,rcons,num)
                    string, lines = PS.sfit_anafit(dfit,key)
                    print_string(string,9)
                    # save lines for KMC
                    LINES.update({(anatype,rctype,case,rcname,direc):line for anatype,line in lines.items()})
                    # save data for plotting
                    if case == 0: d4plot[      rctype] = (khu,dfit)
                    else        : d4plot["anh"+rctype] = (khu,dfit)
                    plotdata.update(manage_data_for_plot_rcons(rcname,direc,ltemp,d4plot,units))
                if plotfile is not None and plotdata != {}: write_plotfile(plotfile,plotdata)

    # create loop list
    loop = [(rctype,case) \
            for case    in [0,1] \
            for rctype  in "tst,cvt,tstzct,cvtzct,tstsct,cvtsct".split(",") \
           ]
    # Print lines for KMC
    string  = ""
    string += "=================\n"
    string += " FITTING SUMMARY \n"
    string += "=================\n"
    string += "\n"
    for rctype,case in loop:
        lines = ""
        for anatype in [1,2,3,4]:
            for rcname in targets:
                for direc in "fw,bw".split(","):
                    key = (anatype,rctype,case,rcname,direc)
                    if key not in LINES.keys(): continue
                    lines += "  "+LINES[key]
            if lines != "": lines += "\n"
        if lines == "": continue
        if case == 0: string += "  HARMONIC-%s\n"%(PS.rckey2str[rctype])
        else        : string += "  ANHARMONIC-%s\n"%(PS.rckey2str[rctype])
        string += "\n"
        string += lines
        string += "\n"
    print_string(string,3)
#===============================================================#

