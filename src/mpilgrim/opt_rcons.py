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
| Module name:  mpilgrim.opt_rcons  |
| Last Update:  2019/04/03 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*
'''

#--------------------------------------------------#
import time
import os
import sys
import numpy   as np
#--------------------------------------------------#
import names   as PN
import rdwr    as RW
import strings as PS
#--------------------------------------------------#
import common.Exceptions as Exc
import common.partfns    as partfns
#--------------------------------------------------#
from   common.fit2anarc  import fit2anarc
from   common.fncs       import prod_list
from   common.fncs       import print_string
import common.physcons as pc
from   common.physcons   import KB, ML, H, AMU
from   common.Logger     import Logger
from   common.Molecule   import Molecule
#--------------------------------------------------#
from   diverse           import ffchecking
from   diverse           import get_input_data
from   diverse           import status_check
from   plotting          import manage_data_for_plot_rcons
from   plotting          import write_plotfile
#--------------------------------------------------#




#==========================================================#
#               Get data from existing files               #
#==========================================================#
def target2list(target):
    # if target is a string --> put it in a list
    if type(target) == type(""): target = [target]
    # (1) target is None??
    if (target is None): return None
    # (2) None in list
    if (None in target): return None
    # (3) no elements
    if len(target)==0  : return None
    # (4) any empty element
    if "" in target    : return None
    # return
    return target
#----------------------------------------------------------#
def data_target(target,dall,ltemp):
    # convert to list
    target = target2list(target)
    # initialize
    tot_V0  = 0.0
    tot_V1  = 0.0
    tot_Q   = [1.0 for T in ltemp]
    tot_ANH = [1.0 for T in ltemp]
    trg_ANH = []
    # everything ok?
    if target is None: return None, None, None, tot_ANH, trg_ANH
    # read data
    for trg in target:
        ctc, itc =  PN.name2data(trg)
        # Get key
        if itc is None: key = PN.struckey(ctc,"msho")
        else          : key = PN.struckey(ctc, itc  )
        # look for data
        try   : V0,V1,pfns = dall["pfn"][key]
        except: raise Exc.NoData
        # add energies to total
        tot_V0 += V0
        tot_V1 += V1
        # update total partition function
        tot_Q   = [q_i*q_j for q_i,q_j in zip(tot_Q,pfns)]
        # anharmonicity
        if ctc in dall["anh"].keys():
           anh     = dall["anh"][ctc]
           tot_ANH = [r1*r2 for r1,r2 in zip(tot_ANH,anh)]
           trg_ANH.append(trg)
    # return data
    return tot_V0, tot_V1, tot_Q, tot_ANH, trg_ANH
#----------------------------------------------------------#
def get_TSratios(ctc,itcs,V1TS,QTS,ltemp,dall):
    '''
    ratio = QTS(i)/QTS * exp(-dE(i) * beta)
    '''
    exception = Exc.LostConformer(Exception)
    tst_ratios = {}
    for itc,weight in itcs:
        key = PN.struckey(ctc,itc)
        if key not in dall["pfn"].keys():
           exception._var = key
           raise exception
        V0_i, V1_i, Q_i = dall["pfn"][key]
        ratio = [0.0 for T in ltemp]
        dE = (V1_i-V1TS)
        for idx,T in enumerate(ltemp):
            ratio[idx] = weight*Q_i[idx]/QTS[idx] * np.exp(-dE/KB/T)
        tst_ratios[itc] = ratio
    # get individual data
    return tst_ratios
#----------------------------------------------------------#
def get_corrfactors(ltemp,ctc,itcs,dall,tst_ratios):
    '''
    k^(x) = 1/(hbeta) * QTS/QR * exp(dE*beta) * GAMMA^(x)
    GAMMA^(x) = sum_i weight_i * QTS_i/QTS * exp(dE_i*beta) * Gamma_i =
              = sum_i ratioTST_i * Gamma_i
    '''
    cf = {}
    cf["cvt"   ] = {}
    cf["cvtzct"] = {}
    cf["cvtsct"] = {}
    cf["tstzct"] = {}
    cf["tstsct"] = {}
    for itc,weight in itcs:
        tsname  = PN.struckey(ctc,itc)
        # individual coefs
        cvt    = dall.get("cvt"   ,{}).get(tsname,None)
        zct    = dall.get("zct"   ,{}).get(tsname,None)
        sct    = dall.get("sct"   ,{}).get(tsname,None)
        cagtst = dall.get("cagtst",{}).get(tsname,None)
        cagcvt = dall.get("cagcvt",{}).get(tsname,None)
        # correction
        try   : cf_cvt    = np.array(prod_list((    cvt,      )))
        except: cf_cvt    = None
        try   : cf_cvtzct = np.array(prod_list((zct,cvt,cagcvt)))
        except: cf_cvtzct = None
        try   : cf_cvtsct = np.array(prod_list((sct,cvt,cagcvt)))
        except: cf_cvtsct = None
        try   : cf_tstzct = np.array(prod_list((zct,    cagtst)))
        except: cf_tstzct = None
        try   : cf_tstsct = np.array(prod_list((sct,    cagtst)))
        except: cf_tstsct = None
        # save
        cf["cvt"   ][itc] = cf_cvt
        cf["cvtzct"][itc] = cf_cvtzct
        cf["cvtsct"][itc] = cf_cvtsct
        cf["tstzct"][itc] = cf_tstzct
        cf["tstsct"][itc] = cf_tstsct
    # get total
    for rc in cf.keys():
        tot_cf = np.array([0.0 for T in ltemp])
        fail   = False
        for itc,weight in itcs:
            tst_ratio = tst_ratios[itc]
            cf_i      = cf[rc][itc]
            if cf_i is None: fail = True; continue
            tot_cf +=  tst_ratio*cf_i
        if fail: cf[rc]["total"] = None
        else   : cf[rc]["total"] = tot_cf
    return cf
#==========================================================#


#==========================================================#
def get_mass_ch(target,dctc,dimasses):
    '''
    target may be a compound or a list of compounds.
    If a list ==> returns the sum
    '''
    # initialize
    tot_mass = 0.0
    tot_ch   = 0
    # target to list
    target = target2list(target)
    if target is None: return None, None
    # get mass and charge!
    for ctc in target:
        ctc = ctc.split(".")[0]
        if ctc in dctc.keys():
           # generate instance of Molecule
           molecule = Molecule()
           # gts files associated to this ctc
           itc,weight = dctc[ctc]._itcs[0]
           gtsfile    = dctc[ctc].gtsfile(itc)
           # isotopic modification
           diso = dctc[ctc]._diso
           if   itc in diso.keys(): imod = diso[itc]
           elif "*" in diso.keys(): imod = diso["*"]
           else                   : imod = None
           # read gts
           molecule.set_from_gts(gtsfile)
           # apply iso
           molecule.apply_imods(imod,dimasses)
           # mass and charge
           tot_mass += float(molecule._mass)
           tot_ch   += int(molecule._ch)
        else:
           tot_mass,tot_ch = None, None
    return tot_mass,tot_ch
#----------------------------------------------------------#
def check_balance(dctc,dimasses,Rs,TS,Ps):

    # get mass and charge
    massRs,chRs = get_mass_ch(Rs,dctc,dimasses)
    massTS,chTS = get_mass_ch(TS,dctc,dimasses)
    massPs,chPs = get_mass_ch(Ps,dctc,dimasses)
 
    # same mass and charge?
    conserved = True
    if None not in [chRs,chTS]: 
       if chRs != chTS                  : conserved = False
       if abs(massRs-massTS) > 0.001/AMU: conserved = False
    if None not in [chPs,chTS]: 
       if chPs != chTS                  : conserved = False
       if abs(massPs-massTS) > 0.001/AMU: conserved = False
    if None not in [chRs,chPs]: 
       if chRs != chPs                  : conserved = False
       if abs(massRs-massPs) > 0.001/AMU: conserved = False

    # generate table
    string  = "Conservation of charge and mass:\n"
    string += "\n"
    string += "   ----------------------------------------\n"
    string += "                     | charge | mass (amu) \n"
    string += "   ----------------------------------------\n"
    if chRs is not None:
       string += "    reactant(s)      |  %3i   | %10.3f \n"%(chRs,massRs*AMU)
    if chTS is not None:
       string += "    transition state |  %3i   | %10.3f \n"%(chTS,massTS*AMU)
    if chPs is not None:
       string += "    product(s)       |  %3i   | %10.3f \n"%(chPs,massPs*AMU)
    string += "   ----------------------------------------\n"
    if not conserved:
       string += "    charge and/or mass NOT conserved!!\n"

    # Return  data
    num_none = [chRs,chTS,chPs].count(None)
    if num_none >= 2: return ""    , True
    else            : return string, conserved
#----------------------------------------------------------#
def basicinfo_reaction(Rs,TS,Ps,dctc,dimasses):
    # Equation of the reaction
    string  = ""
    string += "  equation   : "
    string += " + ".join(Rs)
    string += " --> "
    if TS is not None: string += TS
    string += " --> "
    string += " + ".join(Ps)
    string += "\n\n"
    # Identify reactants, transition state and products
    string += "  reactant(s)     : %s\n"%(" + ".join(Rs))
    if TS is not None: string += "  transition state: %s\n"%TS
    string += "  product(s)      : %s\n"%(" + ".join(Ps))
    string += "\n"
    if len(Rs) == 2 and Rs[0] == Rs[1]:
        string += "  The two reactants are equal.\n"
        string += "  Forward rate constants will be multiplied by 2\n"
        string += "\n"
    if len(Ps) == 2 and Ps[0] == Ps[1]:
        string += "  The two products are equal.\n"
        string += "  Backward rate constants will be multiplied by 2\n"
        string += "\n"
    print_string(string,5)

    # get length of name
    not_in_dctc = []
    data = []
    ml   = 4
    for target in Rs+[TS]+Ps:
        if target is None: continue
        # ctc and itc name
        if "." in target: ctc,itc = target.split(".")
        else            : ctc,itc = target, None
        # list of itcs
        if itc is None:
           if ctc in dctc.keys(): itcs = dctc[ctc]._itcs
           else                 : itcs = None
        else                    : itcs = [(itc,1)]
        # save data?
        if itcs is None:
           not_in_dctc.append(target)
           continue
        data.append( (target,itcs) )
        if len(target) > ml: ml = len(target)

    # anything not defined in dctc?
    if len(not_in_dctc) != 0:
       string = "The following compounds are not defined in '%s':\n"%PN.IFILE1
       for target in not_in_dctc:
           string += "   * %s\n"%target
           if target in Rs: Rs = []
           if target in Ps: Ps = []
       print_string(string,7)

    # head of table and division line
    head = (" %%%is | conformer | weight"%ml)%"Name"
    division = "-"*len(head)
    # listing conformers and weights
    string  = "Conformational flexibility:\n"
    string += "\n"
    string += "   "+division+"\n"
    string += "   "+head+"\n"
    string += "   "+division+"\n"
    for target,itcs in data:
        for idx,(itc,weight) in enumerate(itcs):
            if idx == 0: string += "   "+(" %%%is "%ml)%target+ "|    %3s    |  %2i  \n"%(itc,weight)
            else       : string += "   "+(" %%%is "%ml)%""    + "|    %3s    |  %2i  \n"%(itc,weight)
        string += "   "+division+"\n"
    print_string(string,7)

    # Conservation of mass/ch
    string,conserved = check_balance(dctc,dimasses,Rs,TS,Ps)
    print_string(string,7)

    return conserved, not_in_dctc
#----------------------------------------------------------#
def energies_and_pfns(Rs,TS,Ps,dall,ltemp):
    trganh = []
    # data reactant(s)
    V0R , V1R , QR , ANHR , trganhR  = data_target(Rs,dall,ltemp)
    # data products(s)
    V0P , V1P , QP , ANHP , trganhP  = data_target(Ps,dall,ltemp)
    # data transition state
    V0TS, V1TS, QTS, ANHTS, trganhTS = data_target(TS,dall,ltemp)
    # save targets with anharmonicity
    if trganhR  is not None: trganh += trganhR
    if trganhTS is not None: trganh += trganhTS
    if trganhP  is not None: trganh += trganhP
    # group data and return
    dataR  = (V0R ,V1R ,QR ,ANHR )
    dataTS = (V0TS,V1TS,QTS,ANHTS)
    dataP  = (V0P ,V1P ,QP ,ANHP )
    return dataR, dataTS, dataP, trganh
#----------------------------------------------------------#
def calculate_eqconst(ltemp,nR,nP,dataR,dataP,wfw=1.0,wbw=1.0):
    (V0R ,V1R ,QR) = dataR
    (V0P ,V1P ,QP) = dataP
    # can be calculated?
    if nR == 0 or nP == 0     : return None
    if None in (V1R,V1P,QR,QP): return None
    # calculate it!
    Keq   = partfns.Qs2Kc(ltemp,QR,QP,V1R,V1P)
    # correct Keq
    Keq   = [ K_i*wfw/wbw for K_i in Keq]
   ## get Gibbs
   #Gibbs = partfns.Kc2GFE(ltemp,Keq,dn=nR-nP)
   ## print string
   #print_string(PS.srcons_keq(ltemp,Keq,Gibbs,nR,nP,wfw,wbw),7)
    return Keq
#----------------------------------------------------------#
def calculate_ratecons(ltemp,nR,dataR,dataTS,corrfactors,weight=1.0,key=None):
    rcons = {}
    # TST rate constant
    V0R ,V1R ,QR = dataR
    V0TS,V1TS,QTS= dataTS
    Kpseudo = partfns.Qs2Kc(ltemp,QR,QTS,V1R,V1TS)
    kTST    = partfns.Kc2rate(ltemp,Kpseudo)
    # correct kTST with weights
    kTST    = [k_i*weight for k_i in kTST]
    # save in dictionary
    rcons["tst"] = kTST
    # other rate contants
    for rctype in corrfactors.keys():
        cf = corrfactors[rctype]["total"]
        if cf is None: continue
        rcons[rctype] = kTST*cf
    return rcons
#----------------------------------------------------------#
def fit2analytic(ltemp,rcons,nR):
    # fit to analytic
    dict_fitting = {}
    if len(ltemp) > 1:
       human_units = pc.ML**(nR-1.0) / pc.SECOND
       rctypes = sorted(rcons.keys())
       rctypes.remove("tst")
       rctypes = ["tst"]+rctypes
       for rctype in rctypes:
           k    = [val*human_units for val in rcons[rctype]]
           dfit = fit2anarc(ltemp,k,log=True)
           # save data
           dict_fitting[rctype] = (k,dfit)
    # return data
    return dict_fitting
#----------------------------------------------------------#
def factors_from_TS(TS,V1TS,QTS,ltemp,dctc,dall):
    # list of itcs
    ctc, itc = PN.name2data(TS)
    if itc is None: itcs = dctc[ctc]._itcs
    else          : itcs = [(itc,1)]
    # get TST ratios
    tst_ratios  = get_TSratios(ctc,itcs,V1TS,QTS,ltemp,dall)
    # print TST ratios
    print_string(PS.srcons_tstratios(ltemp,tst_ratios),7)
    # get total correction factors
    corrfactors = get_corrfactors(ltemp,ctc,itcs,dall,tst_ratios)
    # print correction factors
    print_string(PS.string_corrfactors(ltemp,corrfactors),7)
    # print VTST contributions
    print_string(PS.srcons_vtstratios(ltemp,tst_ratios,corrfactors),7)
    return tst_ratios, corrfactors
#----------------------------------------------------------#
def manage_rcons(rcname,ltemp,nR,dataR,dataTS,corrfactors,wfw,direc,ANH4RC,ANH):
    drcons = {}
    if nR == 0: return drcons
    if direc == "fw":
        print "       Calculation of FORWARD  rate constants"
        print
    else:
        print "       Calculation of BACKWARD rate constants"
        print
    # calculate rcons
    idata = (ltemp,nR,dataR,dataTS,corrfactors,wfw,"%s.%s"%(rcname,direc))
    rcons = calculate_ratecons(*idata)
    # anharmonic rate constants
    if ANH:
        # Calculate anharmonic rate constants
        rctypes = rcons.keys()
        for rctype in rctypes:
            RC     = rcons[rctype]
            ANH_RC = [ki*anhi for ki,anhi in zip(RC,ANH4RC)]
            rcons["anh"+rctype] = ANH_RC
    # print rate constants
    print_string(PS.string_rcons(ltemp,rcons,nR,wfw,case=0),11)
    if ANH: print_string(PS.string_rcons(ltemp,rcons,nR,wfw,case=1),11)
    # print gibbs free energies
    print_string(PS.string_gibbsactenergy(ltemp,rcons,nR,case=0),11)
    if ANH: print_string(PS.string_gibbsactenergy(ltemp,rcons,nR,case=1),11)
    # prepare drcons
    drcons = {"%s.%s.%s"%(k,rcname,direc):v for k,v in rcons.items()}
    #--------------#
    # for plotting #
    #--------------#
    # human units: string and value
    hu = pc.ML**(nR-1.0) / pc.SECOND
    if nR != 1: units = "s^{-1} cc^{%i} molecule^{-%i}"%(nR-1,nR-1)
    else      : units = "s^{-1}"
    # generate plotdata
    d4plot = {key:([ki*hu for ki in k],{}) for key,k in rcons.items()}
    plotdata = manage_data_for_plot_rcons(rcname,direc,ltemp,d4plot,units)
    # return data
    return drcons, plotdata
#----------------------------------------------------------#




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
    unknown = sorted([target for target in targets if target not in dchem.keys()])
    targets = sorted([target for target in targets if target in dchem.keys()])

    if unknown != []:
       print "   The following reactions are not defined:"
       for unkn in unknown: print "      * %s"%unkn
       print

    # read dof
    dall = RW.read_alldata(dof,ltemp)[0]

    #--------------------------------#
    # Calculations for each reaction #
    #--------------------------------#
    for rcname in targets:
        # pof file
        pof = PN.get_pof(dlevel,"rcons",rcname)
        sys.stdout = Logger(pof,"w",True)
        # print title
        print_string(PS.title_rcons(rcname,pof),3)
        # get data
        Rs, TS, Ps = dchem[rcname]
        # Rs = Ps?
        wfw, wbw = 1, 1
        if len(Rs) == 2 and Rs[0] == Rs[1]: wfw *= 2
        if len(Ps) == 2 and Ps[0] == Ps[1]: wbw *= 2
        # basic info
        allright, not_in_dctc = basicinfo_reaction(Rs,TS,Ps,dctc,dimasses)
        if not allright: continue
        # clean data
        if len( set(not_in_dctc).intersection(Rs) ) != 0: Rs = []
        if len( set(not_in_dctc).intersection(Ps) ) != 0: Ps = []
        if TS in not_in_dctc: TS = None
        nR = len(Rs)
        nP = len(Ps)
        if nR == 0:
           print "     Reactant(s) not defined or unknown!"
           print
           continue
        # energies and partition functions
        dataR,dataTS,dataP,trganh = energies_and_pfns(Rs,TS,Ps,dall,ltemp)
        ANH = (len(trganh) != 0)
        (V0R ,V1R ,QR ,ANHR ) = dataR
        (V0TS,V1TS,QTS,ANHTS) = dataTS
        (V0P ,V1P ,QP ,ANHP ) = dataP
        dataR  = (V0R ,V1R ,QR)
        dataTS = (V0TS,V1TS,QTS)
        dataP  = (V0P ,V1P ,QP)
        # anharmonicity??
        if ANH:
           # anharmonicity factors
           ANH4KEQ = [ anhP /anhR for anhR,anhP  in zip(ANHR,ANHP )]
           ANH4RCF = [ anhTS/anhR for anhR,anhTS in zip(ANHR,ANHTS)]
           ANH4RCB = [ anhTS/anhP for anhP,anhTS in zip(ANHP,ANHTS)]
           bools = [nR*nP != 0,nR!=0,nP!=0]
           print_string(PS.srcons_anhtable(trganh,ltemp,ANH4KEQ,ANH4RCF,ANH4RCB,bools),7)
        else:
           ANH4KEQ = None
           ANH4RCF = None
           ANH4RCB = None
        # print energy table
        print_string(PS.srcons_relenergies(dctc,dall,Rs,TS,Ps,V0R,V1R),7)
        #-----------------------#
        # equilibrium constants #
        #-----------------------#
        if nR*nP != 0:
           # harmonic
           Keq  = calculate_eqconst(ltemp,nR,nP,dataR,dataP,wfw,wbw)
           print_string(PS.srcons_keq(ltemp,Keq,nR,nP,wfw,wbw,0),7)
           # anharmonic
           if ANH:
              aKeq = [Ki*anhi for Ki,anhi in zip(Keq,ANH4KEQ)]
              print_string(PS.srcons_keq(ltemp,aKeq,nR,nP,wfw,wbw,1),7)
        #-----------------------#
        #    rate  constants    #
        #-----------------------#
        if TS is not None:
            # TST ratios and correction factors to TST
           idata = (TS,V1TS,QTS,ltemp,dctc,dall)
           tst_ratios, corrfactors = factors_from_TS(*idata)
           # FORWARD and BACKWARD Rate constants
           drcons,plotdata = {},{}
           for direc,num,weight,data,ANH4RC in [("fw",nR,wfw,dataR,ANH4RCF),("bw",nP,wbw,dataP,ANH4RCB)]:
               if num == 0: continue
               idata = (rcname,ltemp,num,data,dataTS,corrfactors,weight,direc,ANH4RC,ANH)
               o1,o2 = manage_rcons(*idata)
               drcons.update(o1)
               plotdata.update(o2)
           # save data for plotting
           if plotfile is not None and plotdata != {}: write_plotfile(plotfile,plotdata)
        else:
           drcons = {}
        #-------------#
        # update dall #
        #-------------#
        print "   Updating data file: %s"%dof
        dall = RW.read_alldata(dof,ltemp)[0]
        if "rcons" not in dall.keys(): dall["rcons"] = {}
        dall["rcons"].update(drcons)
        RW.write_alldata(dof,ltemp,dall)
        print
        #-----------------------#
        # Reactants == Products #
        #-----------------------#
        if sorted(Rs) == sorted(Ps):
           print "  =============================================="
           print "  ||  Reactants and products are identical.   ||"
           print "  ||  Rate constants will be multiplied by 2. ||"
           print "  =============================================="
           print
           rcons = {k.split(".")[0]:v for k,v in drcons.items() if k.endswith("fw")}
           for key in rcons.keys():
               rcons[key] = np.array(rcons[key]) * 2
           print_string(PS.string_rcons(ltemp,rcons,nR,case=0),11)
           if ANH: print_string(PS.string_rcons(ltemp,rcons,nR,case=1),11)
#===============================================================#

