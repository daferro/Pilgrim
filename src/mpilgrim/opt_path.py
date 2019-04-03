#!/usr/bin/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.opt_path   |
| Last Update:  2019/03/29 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*
'''

#--------------------------------------------------#
import copy
import gc
import os
import readline
import sys
import time
import numpy             as np
import shutil
#--------------------------------------------------#
import names      as PN
import strings    as PS
import rdwr       as RW
#--------------------------------------------------#
from   itf        import get_spc_fnc
from   diverse    import get_input_data
from   diverse    import status_check
from   diverse    import ffchecking
from   diverse    import find_label_in_rst
#--------------------------------------------------#
from   plotting   import manage_data_for_plot_mep
from   plotting   import manage_data_for_plot_cvt
from   plotting   import manage_data_for_plot_sct
from   plotting   import write_plotfile
from   exceptions import deal_with_exception
#--------------------------------------------------#
import common.fncs       as fncs
import common.adipot     as ap
import common.cvtst      as cv
import common.files      as ff
import common.internal   as intl
import common.physcons   as pc
import common.cag        as cag
import common.sct        as sct
import common.steepdesc  as sd
import common.Exceptions as Exc
import common.ispe       as ispe
#--------------------------------------------------#
from   common.Logger     import Logger
from   common.Molecule   import Molecule
from   common.steepdesc  import TSLABEL
from   common.criteria   import EPS_MEPS
#--------------------------------------------------#



WARNINGS = []

#===============================================================#
def get_itargets(targets,dpath,dctc):
    '''
    get individual targets
    '''
    # Targets?
    if (len(targets) == 0) or ("*" in targets): targets = dpath.keys()

    # generate list of individual targets
    itargets = []
    for target in targets:
        ctc, itc = PN.name2data(target)
        # check ctc
        if ctc not in dctc.keys() :
           print "   * '%s' not in '%s'"%(target,PN.IFILE1)
           print
           continue
        if ctc not in dpath.keys():
           print "   * '%s' not in '%s'"%(target,PN.IFILE3)
           print
           continue
        # list of itcs for the ctc
        itclist  = [itc_i for itc_i,weight_i in dctc[ctc]._itcs]
        # add itc
        if   itc is None    : itargets += [PN.struckey(ctc,itc) for itc,weight in dctc[ctc]._itcs]
        elif itc in itclist : itargets += [PN.struckey(ctc,itc)]
    return sorted(itargets)
#---------------------------------------------------------------#
def in_interval(s,si,sf,eps=1e-8):
    return si-eps <= s <= sf+eps
#---------------------------------------------------------------#
def get_masses(target,dctc,dimasses):
    ctc, itc = PN.name2data(target)
    diso     = dctc[ctc]._diso
    if   itc in diso.keys(): imod = diso[itc]
    elif "*" in diso.keys(): imod = diso["*"]
    else                   : imod = None
    if imod is None: return None
    gtsTS  = dctc[ctc].gtsfile(itc)
    TS     = Molecule()
    TS.set_from_gts(gtsTS)
    TS.apply_imods(imod,dimasses)
    masses = list(TS._masses)
    return masses
#===============================================================#


#===============================================================#
def compare_tpath(tpath,tpath2,rstfile,eps=1e-8):
    exception = Exc.RstDiffVar(Exception)
    exception._rst = rstfile
    if tpath2 is None: return
    the_vars = ("method","mu","ds","hsteps","cubic")
    for idx,var in enumerate(the_vars):
        exception._var = var
        var1, var2 = tpath[idx],tpath2[idx]
        if type(var1) != type(var2): raise exception
        if type(var1) == type(""):
            if var1 != var2: raise exception
        if type(var1) == type(1.0):
            if abs(var1-var2) > eps: raise exception
#---------------------------------------------------------------#
def compare_tcommon(tcommon,tcommon2,rstfile,eps=1e-8):
    exception = Exc.RstDiffVar(Exception)
    exception._rst = rstfile
    if tcommon2 is None: return
    the_vars = ("ch","mtp","atonums","masses","mu")
    for idx,var in enumerate(the_vars):
        exception._var = var
        var1, var2 = tcommon[idx],tcommon2[idx]
        if type(var1) == type(1) or type(var1) == type(1.0):
            var1 = [var1]
            var2 = [var2]
        for ii,jj in zip(var1,var2):
            if abs(ii-jj) > eps: raise exception
#===============================================================#

#===============================================================#
def calc_mep(itarget,gtsTS,pathvars,tsoftware,TMP,decrease=False):
    '''
     if decrease = True, MEP is reduced hsteps steps
    '''
    # data in name
    ctc, itc = PN.name2data(itarget)
    # calculate path
    tcommon,drst,pathvars = obtain_mep(itarget,gtsTS,pathvars,tsoftware,TMP)
    # decrease mep??
    if decrease:
       sbw1,sfw1 = pathvars._sbw,pathvars._sfw
       pathvars.decrease_svals()
       sbw2,sfw2 = pathvars._sbw,pathvars._sfw
       print "    MEP will be reduced to check SCT convergence:"
       print
       print "        sbw: %+8.4f --> %+8.4f bohr"%(sbw1,sbw2)
       print "        sfw: %+8.4f --> %+8.4f bohr"%(sfw1,sfw2)
       print
       drst = {label:data for label,data in drst.items() if in_interval(data[0],sbw2,sfw2)}
    # DLEVEL??
    if pathvars._dlevel is not None:
       print "    Applying Dual-Level..."
       print
       dlevel_xy = [(find_label_in_rst(x,drst)[0],y) for x,y in pathvars._dlevel.items()]
       dlevel_xy = [(x,y) for x,y in dlevel_xy if x is not None]
       # Print points (sorted by s value)
       dummy = [(drst[xx][0],idx) for idx,(xx,yy) in enumerate(dlevel_xy)]
       dummy.sort()
       for s_i,idx_i in dummy:
           xx_i, yy_i = dlevel_xy[idx_i]
           print "             %+8.4f bohr (%-6s) --> %.7f hartree"%(s_i,xx_i,yy_i)
       print
       # interpolation
       drst, points, xx, yyll, yyhl = ispe.ispe(tcommon,drst,dlevel_xy,tension=0.0)
       tdleveldata = (points,xx,yyll,yyhl)
       # table new values
       fncs.print_string(PS.smep_tableDLEVEL(drst,tdleveldata,pathvars._eref),4)
    else:
       # Print table
       fncs.print_string(PS.smep_table(drst,pathvars._eref),4)
    return tcommon, drst, pathvars
#---------------------------------------------------------------#
def onesidemep(ivars,rstfile,drst):
    points = []
    for geom in sd.steepest(*ivars):
        (point,si), E, xms, gms, Fms, t = geom
        points.append(point)
        if point not in drst.keys():
           drst[point] = (si,E,xms,gms,Fms,None,None,t)
           ff.write_rst_add(rstfile,point,drst[point])
    # check convergence of last point
    if len(points) > 2:
       masses = ivars[3]
       mu     = ivars[4][1]
       epse   = ivars[4][5]
       epsg   = ivars[4][6]
       sA, EA, xA, gA = drst[points[-2]][0:4]
       sB, EB, xB, gB = drst[points[-1]][0:4]
       # check convergence of MEP
       dE   = abs(EB-EA)
       ngB  = fncs.norm(fncs.ms2cc_g(gB,masses,mu))
       conv = (dE  < epse) and (ngB < epsg)
       s_y = "%+6.3f"%drst[points[-2]][0]
       s_z = "%+6.3f"%drst[points[-1]][0]
       tt1 = "|E(s=%s)-E(s=%s)|"%(s_z,s_y)
       tt2 = "|grad(s=%s)|"%(s_z)
       ml = max(max([len(tt1),len(tt2)]),20)
       line1 = ("%%-%is = %%.2E hartree     "%ml)%(tt1,dE)
       line2 = ("%%-%is = %%.2E hartree/bohr"%ml)%(tt2,ngB)
       if (dE  < epse): line1 += " < %.2E  YES"%epse
       else           : line1 += " > %.2E  NO"%epse
       if (ngB < epsg): line2 += " < %.2E  YES"%epsg
       else           : line2 += " > %.2E  NO"%epsg
    else:
       conv  = False
       line1 = ""
       line2 = ""
    tconv = (conv,line1,line2)
    return (drst,points,tconv)
#---------------------------------------------------------------#
def obtain_mep(target,gtsTS,pathvars,tsoftware,TMP):
    ctc, itc = PN.name2data(target)

    # Files
    rstfile = PN.get_rst(ctc,itc)
    xyzfile = PN.get_rstxyz(ctc,itc)

    # print
    software, tes = tsoftware
    fncs.print_string(PS.smep_init(target,software,pathvars._paral,pathvars.tuple_first(),\
                 pathvars.tuple_sdbw(),pathvars.tuple_sdfw()),4)
    fncs.print_string(PS.smep_ff(TMP,PN.DIR4,PN.DIR5,rstfile,xyzfile),4)

    # data for single-point calculation
    spc_args  = (tes,TMP,False)
    spc_fnc   = get_spc_fnc(software)

    # read rst
    tpath2, tcommon2, drst = ff.read_rst(rstfile)
    fncs.print_string(PS.smep_rst(rstfile,drst),4)
    
    # correct MEP direction?
    if TSLABEL in drst.keys():
       ii_s, ii_V, ii_x, ii_g, ii_F, ii_v0, ii_v1, ii_t = drst[TSLABEL]
       ii_ic, ii_sign = pathvars._fwdir
       if not intl.ics_correctdir(ii_x,ii_v0,ii_ic,ii_sign,tcommon2[3],tcommon2[4]):
           print "    'fwdir' variable differs from MEP direction in rst file!"
           print "        * modifying rst internal dictionary..."
           new_drst = {}
           for key in drst.keys():
               ii_s, ii_V, ii_x, ii_g, ii_F, ii_v0, ii_v1, ii_t = drst[key]
               ii_s = -ii_s
               if ii_v0 is not None: ii_v0 = [-ii for ii in ii_v0]
               if ii_v1 is not None: ii_v1 = [-ii for ii in ii_v1]
               if "bw" in key: newkey = key.replace("bw","fw")
               else          : newkey = key.replace("fw","bw")
               new_drst[newkey] = (ii_s, ii_V, ii_x, ii_g, ii_F, ii_v0, ii_v1, ii_t)
           drst = new_drst
           del new_drst
           print "        * rewriting rst file..."
           ff.write_rst(rstfile,tpath2,tcommon2,drst)

    # Extension of MEP in rst is bigger
    if drst != {}:
       lbw,lfw,sbw,sfw,Ebw,Efw = sd.rstlimits(drst)
       pathvars._sbw = min(pathvars._sbw,sbw)
       pathvars._sfw = max(pathvars._sfw,sfw)

    # Read gts
    ts = Molecule()
    ts.set_from_gts(gtsTS)
    # scaling of frequencies
    ts.set_fscal(pathvars._freqscal)
    # apply iso mod
    if pathvars._masses is not None: ts.mod_masses(pathvars._masses)
    # setup
    ts.setup(mu=pathvars._mu)
    ts.ana_freqs(case="cc")

    fncs.print_string(PS.smep_ts(ts),4)

    tcommon = (ts._ch,ts._mtp,ts._atnums,ts._masses,ts._mu)
    compare_tpath(pathvars.tuple_rst(),tpath2,rstfile)
    compare_tcommon(tcommon,tcommon2,rstfile)

    #------------#
    # First step #
    #------------#
    print
    print "    Performing first step of MEP..."
    print
    inputvars = (ts._xcc,ts._gcc,ts._Fcc,ts._symbols,ts._masses,pathvars.tuple_first(),\
                 spc_fnc,spc_args,drst,pathvars._paral)
    (xms,gms,Fms),(v0,v1),(xms_bw,xms_fw) = sd.mep_first(*inputvars)
    s1bw = -float(pathvars._ds)
    s1fw = +float(pathvars._ds)
    fncs.print_string(PS.smep_first(ts._symbols,ts._xms,v0,v1),8)

    # write rst file
    if TSLABEL not in drst.keys():
        drst[TSLABEL] = (0.0,ts._V0,xms,gms,Fms,v0,v1,None)
        ff.write_rst_head(rstfile,pathvars.tuple_rst(),tcommon)
        ff.write_rst_add(rstfile,TSLABEL,drst[TSLABEL])

    #------------#
    # The MEP    #
    #------------#
    print
    print "    Calculating MEP..."
    print
    print "       * REMEMBER: data of each step will be saved at %s"%rstfile
    print "                   a summary will be printed when finished"
    # preparation
    xcc_bw = fncs.ms2cc_x(xms_bw,ts._masses,pathvars._mu)
    xcc_fw = fncs.ms2cc_x(xms_fw,ts._masses,pathvars._mu)
    args_bw = ((xcc_bw,s1bw,ts._symbols,ts._masses,pathvars.tuple_sdbw(),\
                spc_fnc,spc_args,drst,ts._Fms,"bw%i") , rstfile,drst)
    args_fw = ((xcc_fw,s1fw,ts._symbols,ts._masses,pathvars.tuple_sdfw(),\
                spc_fnc,spc_args,drst,ts._Fms,"fw%i") , rstfile,drst)
    # execution
    if pathvars._paral:
        import multiprocessing
        pool = multiprocessing.Pool()
        out  = [pool.apply_async(onesidemep,args=args) for args in [args_bw,args_fw]]
        drstbw,pointsbw,convbw = out[0].get()
        drstfw,pointsfw,convfw = out[1].get()
        del out
        # clean up pool
        pool.close()
        pool.join()
    else:
        drstbw,pointsbw,convbw = onesidemep(*args_bw)
        drstfw,pointsfw,convfw = onesidemep(*args_fw)
    # update drst
    print "       * FINISHED!"
    print
    drst.update( drstbw )
    drst.update( drstfw )
    points = [TSLABEL]+pointsbw+pointsfw
    # empty variables
    del drstbw,pointsbw
    del drstfw,pointsfw
    # Rewrite rst
    print "       * (re)writing file: %s (sorted version)"%rstfile
    ff.write_rst(rstfile,pathvars.tuple_rst(),tcommon,drst)
    print
   ## restrict drst to points
   #if restrict: drst = {point:drst[point] for point in points}
    # Get limits of rst
    lbw,lfw,sbw,sfw,Ebw,Efw = sd.rstlimits(drst)
    convbw,l1bw,l2bw = convbw
    convfw,l1fw,l2fw = convfw
    if l1bw+l1fw != "":
        print "       * EPS criteria (epse and epsg):"
        print
        if l1bw != "": print "         %s"%l1bw
        if l2bw != "": print "         %s"%l2bw
        if l1fw != "": print "         %s"%l1fw
        if l2fw != "": print "         %s"%l2fw
        print
        if convbw:
           pathvars.converged_in_bw(sbw)
           print "         CRITERIA FULFILLED in backward dir.!"
           print "         path stopped at sbw = %+8.4f bohr"%sbw
           print
        if convfw:
           pathvars.converged_in_fw(sfw)
           print "         CRITERIA FULFILLED in forward dir.!"
           print "         path stopped at sfw = %+8.4f bohr"%sfw
           print
    # write molden file
    print "       * writing file: %s"%xyzfile
    ff.rst2xyz(rstfile,xyzfile,onlyhess=True)
    print 
    # reference energy
    if pathvars._eref is None: pathvars._eref = Ebw
    # return data
    return tcommon, drst, pathvars
#---------------------------------------------------------------#
def obtain_drp():
    pass
#===============================================================#


#===============================================================#
def calc_coefs(itarget,tcommon,drst,pathvars,ltemp,plotfile=None):
    # initialize dcfs
    dcfs = {}
    # sorted points
    points = sd.sorted_points(drst,hess=True)
    # Calculate adiabatic potential
    dMols, Vadi, pathvars = obtain_adipot(tcommon,drst,pathvars)
    # Calculate CVT correction factor
    if pathvars._cvt == "yes":
       dcfs, lscvt, gibbs,lnew = obtain_cvt(dMols,points,Vadi,ltemp,pathvars,dcfs=dcfs)
    else:
       lscvt = None
    # Calculate CAG correction factor
    dcfs = obtain_cag(Vadi,ltemp,lscvt,dcfs=dcfs)
    # Calculate SCT correction factor
    if pathvars._sct == "yes":
       # v1 vector along MEP
       if   pathvars._v1mode == "grad": dv1 = sct.get_numv1(drst)
       elif pathvars._v1mode == "hess": dv1 = {}
       # calculate SCT coefficient
       dcfs, tplot_sct = obtain_sct(dMols,points,Vadi,ltemp,dv1,pathvars,dcfs=dcfs)
    # save data for plotting
    if plotfile is not None:
       plotdata = {}
       plotdata.update(manage_data_for_plot_mep(itarget,drst,pathvars._eref,Vadi))
       if "cvt" in dcfs.keys():
         plotdata.update(manage_data_for_plot_cvt(itarget,ltemp,Vadi.xx(),(gibbs,lnew),lscvt,dcfs["cvt"]))
       if "sct" in dcfs.keys():
         plotdata.update(manage_data_for_plot_sct(itarget,dcfs["zct"],dcfs["sct"],*tplot_sct))
       write_plotfile(plotfile,plotdata)
    # Delete variables
    del dMols
    del Vadi
    # Return data
    return dcfs, pathvars
#---------------------------------------------------------------#
def obtain_adipot(tcommon,drst,pathvars):
    print
    print "    Calculating adiabatic potential..."
    print 
    ics = pathvars.get_ics()
    if pathvars._useics == "yes" and (ics is None or len(ics) == 0):
       raise Exc.NoICS(Exception)
    # calculate adiabatic potential
    idata = (tcommon,drst,pathvars._eref,ics,pathvars._useics,pathvars._lowfq,pathvars._freqscal)
    dMols, ccVadi, icVadi, tuple_cc, tuple_ic, listV0 = ap.path2vadi(*idata)
    # expand tuples
    data_x, lcc_frqs , lcc_tzpe = tuple_cc
    data_x, lic_frqs , lic_tzpe = tuple_ic
    # check data
    fncs.print_string(PS.sadipot_ics(ics,pathvars._useics),8)
    if pathvars._useics == "yes":
       ok1 = ap.ccvsic_checknfreqs(lcc_frqs,lic_frqs)
       ok2 = ap.ccvsic_checkts(data_x,lcc_frqs,lic_frqs)
       ok3 = ap.ccvsic_checkzpe(lcc_tzpe,lic_tzpe)
       checks = (ok1,ok2,ok3)
       fncs.print_string(PS.sadipot_checks(checks),8)
       # select spline
       if not ok1 or not ok2          : Vadi, pathvars._useics = ccVadi, "no"
       elif pathvars._useics == "yes" : Vadi, pathvars._useics = icVadi, "yes"
       else                           : Vadi, pathvars._useics = ccVadi, "no"
    else:
       Vadi, pathvars._useics = ccVadi, "no"
    # setup spline
    Vadi.setup()
    data4print = (Vadi.xx(), Vadi.yy(), Vadi._sAG, Vadi._VAG, listV0, lcc_tzpe, lic_tzpe, pathvars._eref)
    fncs.print_string(PS.sadipot_table(*data4print),8)
    # print freqs
    print "      Vibrational frequencies summary: cc (ic) [in cm^-1]:"
    fncs.print_string(PS.sadipot_freqs(Vadi.xx(),lcc_frqs,lic_frqs),8)
    return dMols, Vadi, pathvars
#---------------------------------------------------------------#
def obtain_cvt(dMols,points,VadiSpl,temps,pathvars,si=-float("inf"),sj=+float("inf"),dcfs={}):
    print
    print "    Calculating CVT variational coefficient..."
    print 
    useics = pathvars._useics
    if len(temps) == 0: raise Exc.NoTemps(Exception)
    # Only points between si and sj
    points = [pp for pp in points if si <= dMols[pp][0] <= sj]
    lcvt_s, lcvt_gamma, gibbs_matrix, gibbsTS, lnew = cv.get_cvt(dMols,points,VadiSpl,temps,useics)
    # print gibbs
    svals = [dMols[point][0] for point in points]
    fncs.print_string(PS.scvt_gibbs(svals,temps,gibbs_matrix.copy(),pathvars,gibbsTS),8)
    # print cvt coefs
    fncs.print_string(PS.scvt_coefs(lcvt_s, lcvt_gamma, temps),8)
    # save data
    dcfs["cvt"] = lcvt_gamma
    return dcfs, lcvt_s, gibbs_matrix, lnew
#---------------------------------------------------------------#
def obtain_cag(Vadi,ltemp,lscvt=None,dcfs={}):
    print
    print "    Calculating CAG coefficient..."
    print 
    dE_cagtst, cagtst = cag.calc_cag(ltemp,Vadi)
    dcfs["cagtst"] = cagtst
    if lscvt is not None:
       dE_cagcvt, cagcvt = cag.calc_cag(ltemp,Vadi,lscvt)
       dcfs["cagcvt"] = cagcvt
    else:
       dE_cagcvt = None
       cagcvt    = None
    fncs.print_string(PS.scag_table(ltemp,dE_cagtst,cagtst,dE_cagcvt,cagcvt),8)
    return dcfs
#---------------------------------------------------------------#
def obtain_sct(dMols,points,VadiSpl,temps,dv1,pathvars,dcfs={}):
    print
    print "    Calculating SCT transmission coefficient..."
    print 
    # data from  pathvars
    useics = pathvars._useics
    v1mode = pathvars._v1mode
    # E0 value
    if pathvars._e0 is None:
       V1bw = pathvars._eref + VadiSpl.get_alpha()[1]
       V1fw = pathvars._eref + VadiSpl.get_omega()[1]
       if pathvars._V1R is not None: E0bw = pathvars._V1R
       else                        : E0bw = V1bw
       if pathvars._V1P is not None: E0fw = pathvars._V1P
       else                        : E0fw = V1fw
       E0 = max(E0bw,E0fw) - pathvars._eref
    else:
       E0 = pathvars._e0 - pathvars._eref
    # some checks
    if len(temps) == 0: raise Exc.NoTemps(Exception)
    if useics in ["yes",True]: case = "ic"
    else                     : case = "cc"
    # Part I - Get E0 and VAG
    E0      = sct.get_sct_part1(points,VadiSpl,E0)
    sAG,VAG = VadiSpl.get_max()
    fncs.print_string(PS.ssct_init(E0,VadiSpl,pathvars,v1mode),8)
    # Part II - Calculate tbar, bmfs and mueff
    tuple_part2 = (dMols,points,dv1,case,pathvars._muintrpl)
    svals, lkappa, ltbar, ldtbar, mu, lmueff, toignore = sct.get_sct_part2(*tuple_part2)
    fncs.print_string(PS.ssct_mueff(svals,VadiSpl,lkappa,ltbar,lmueff,toignore),8)
    # Part III - Quantum reaction coordinate
    fncs.print_string(PS.ssct_E0VAG(E0,VAG),8)
    if pathvars._qrc is not None:
       afreq   = pathvars._qrcafreq
       lEquant = pathvars._qrclE
       print "        Quantum reaction coordinate keyword (qrc) activated!"
       print
       mode  = pathvars._qrc[0]+1
       numst = pathvars._qrc[1]
       print "           * reactant mode     : %i (%.2f cm^-1)"%(mode,fncs.afreq2cm(afreq))
       print "           * number of states  : %i"%numst
       print "           * contribution to kappa_SCT from E0 to VAG will"
       print "             be obtained from discrete set of energies"
       print
       print "           * calculating transmission probabilities..."
       print
       qrc_ZCT = sct.get_sct_part3(svals, mu   ,VadiSpl,afreq,lEquant,E0,VAG,temps)
       qrc_SCT = sct.get_sct_part3(svals,lmueff,VadiSpl,afreq,lEquant,E0,VAG,temps)
       nE = len(qrc_SCT[1])
       fncs.print_string(PS.ssct_probs(qrc_SCT[1],qrc_ZCT[2],qrc_SCT[2],qrc_SCT[3]),12)
       print "           * number of included states  : %i"%nE
       print
       kappaI1_zct = qrc_ZCT[0]
       kappaI1_sct = qrc_SCT[0]
    else:
       kappaI1_zct = None
       kappaI1_sct = None
    # Part IV - calculate thetas and probs
    print "        Transmission probabilities for kappa_SCT calculation:"
    print
    weights_ZCT,lE_ZCT,probs_ZCT,rpoints_ZCT,diffs_ZCT = sct.get_sct_part4(svals,mu    ,VadiSpl,E0)
    weights_SCT,lE_SCT,probs_SCT,rpoints_SCT,diffs_SCT = sct.get_sct_part4(svals,lmueff,VadiSpl,E0)
    fncs.print_string(PS.ssct_probs(lE_SCT,probs_ZCT,probs_SCT,rpoints_SCT),8)
    fncs.print_string(PS.ssct_diffs(lE_SCT,diffs_SCT),8)
    # Part V - calculate coefficients
    ZCTdata = sct.get_sct_part5(lE_ZCT,probs_ZCT,weights_ZCT,E0,VAG,temps,kappaI1_zct)
    SCTdata = sct.get_sct_part5(lE_SCT,probs_SCT,weights_SCT,E0,VAG,temps,kappaI1_sct)
    ZCT,lIi_ZCT, RTE_ZCT, INTG_ZCT = ZCTdata
    SCT,lIi_SCT, RTE_SCT, INTG_SCT = SCTdata
    fncs.print_string(PS.ssct_kappa(temps,ZCT,lIi_ZCT,RTE_ZCT,E0,case="zct"),8)
    fncs.print_string(PS.ssct_kappa(temps,SCT,lIi_SCT,RTE_SCT,E0,case="sct"),8)
    #fncs.print_string(PS.ssct_kappa(temps,ZCT,SCT,RTE_ZCT,RTE_SCT,E0),8)
    # save data
    dcfs["zct"] = ZCT
    dcfs["sct"] = SCT
    # data for the plot
    forplot = (svals,lmueff,temps,INTG_ZCT,INTG_SCT,RTE_ZCT,RTE_SCT,E0,VAG)
    return dcfs, forplot
#===============================================================#
def get_path_sctconv(itarget,gtsTS,pathvars,tsoftware,ltemp,TMP,plotfile):
    # assert no DLEVEL will be done (just in case)
    pathvars._dlevel = None
    # go ahead
    convlist_sct = []
    convlist_lim = []
    for step in range(pathvars._sctmns+2):
        if step == 0: decrease = True
        else        : decrease = False
        tcommon,drst,pathvars = calc_mep(itarget,gtsTS,pathvars,tsoftware,TMP,decrease)
        dcfs,pathvars         = calc_coefs(itarget,tcommon,drst,pathvars,ltemp,plotfile)
        # expand data
        lbw, lfw, sbw, sfw, Ebw, Efw = sd.rstlimits(drst)
        SCT = dcfs["sct"]
        del drst
        # Save data only for lower temperature
        convlist_sct.append(SCT[0])
        convlist_lim.append((sbw,sfw))
        # print convergence table
        fncs.print_string(PS.ssct_convergence(convlist_sct,convlist_lim,pathvars._scterr),4)
        # Check convergence
        if len(convlist_sct) > 1:
           SCT_a  = convlist_sct[-1]
           SCT_b  = convlist_sct[-2]
           dif100 = 100*abs(SCT_a-SCT_b)/SCT_a
           if dif100 < pathvars._scterr: break
        # last  step and not converged
        if step+1 == pathvars._sctmns+2:
           print "    WARNING: SCT has not converged but 'sctmns' was reached!"
           print
           break
        # increase variables sbw and sfw
        sbw1,sfw1 = pathvars._sbw, pathvars._sfw
        pathvars.increase_svals(Ebw,Efw)
        sbw2,sfw2 = pathvars._sbw, pathvars._sfw
        if abs(sbw1-sbw2)<EPS_MEPS and abs(sfw1-sfw2)<EPS_MEPS:
           print "    WARNING: MEP cannot be increased and SCT is not converged!"
           print "    Maybe epse and epsg criteria should be modified..."
           print
           break
        else:
           print "    MEP will be increased for SCT convergence:"
           print
           print "        sbw: %+8.4f --> %+8.4f bohr"%(sbw1,sbw2)
           print "        sfw: %+8.4f --> %+8.4f bohr"%(sfw1,sfw2)
           print
           print "   =================================================="
           print
    return dcfs
#---------------------------------------------------------------#
def deal_with_path(target,dlevel,software,ltemp,dctc,pathvars,dtes,dchem,dhighlvl,dimasses):
    dof      = PN.get_dof(dlevel)
    plotfile = PN.get_plf(dlevel)
    # gts file for this TS
    ctc, itc = PN.name2data(target)
    gtsTS    = dctc[ctc].gtsfile(itc)
    # temporal folder
    TMP = PN.TMPi%(target)
    # if exists,remove content
    # (to avoid reading old files from different levels of calculation)
    if os.path.exists(TMP): shutil.rmtree(TMP,ignore_errors=True)
    # split target
    ctc, itc = PN.name2data(target)
    # name of rst file
    rstfile = PN.get_rst(ctc,itc)
    # tuple software
    tes    = dtes.get(software,{}).get(ctc,None)
    tsoftw = (software,tes)
    # internal coordinates
    dics = dctc[ctc]._dics
    if   itc in dics.keys(): ics = dics[itc]
    elif "*" in dics.keys(): ics = dics["*"]
    else                   : ics = None
    # path variables
    pathvars.set_ics(ics) # before setup3!!
    pathvars.apply_specific(itc)
    pathvars.setup1()
    pathvars.setup2()
    pathvars.setup3()
    # Set Eref (from reaction)
    pathvars.set_eref_from_reaction(target,dchem,dof)
    # Quantum reaction coordinate qrc
    pathvars.prepare_qrc(dchem,dctc,dimasses)
    # frequency scaling factor
    pathvars._freqscal = float(dctc[ctc]._fscal)
    # if dlevel --> no convergence and dlevel data
    if dlevel:
       exception = Exc.NoDLEVELdata(Exception)
       pathvars._sctmns = 0
       keydhl = "%s.%s.path"%(ctc,itc)
       if keydhl not in dhighlvl.keys():
          # maybe only TS
          keydhl = "%s.%s"%(dctc[ctc]._root,itc)
          if keydhl in dhighlvl.keys():
             dictE  = {0.0:dhighlvl[keydhl]}
          else:
            global WARNINGS
            WARNINGS.append("No high-level data for %s"%target)
            raise exception
       else: dictE = dhighlvl[keydhl]
       pathvars._dlevel = dictE
    # LOGGER
    pof = PN.get_pof(dlevel,"path",target)
    sys.stdout = Logger(pof,"w",True)
    #string
    fncs.print_string(PS.smep_title(target,pathvars,pof),2)
    #----------------#
    # calculate path #
    #----------------#
    # 1. MEP
    if not pathvars._beyondmep:
        common,drst,pathvars = calc_mep(target,gtsTS,pathvars,tsoftw,TMP,decrease=False)
        dcoefs = {}
        del drst
        # raise error
        raise Exc.OnlyMEP(Exception)
    # 2. MEP + coefs
    else:
       # 2.1 MEP expanded till SCT convergence
       if pathvars.sct_convergence():
          dcoefs = get_path_sctconv(target,gtsTS,pathvars,tsoftw,ltemp,TMP,plotfile)
       # 2.2 Coefs with the current MEP extension
       else:
          tcommon,drst,pathvars = calc_mep(target,gtsTS,pathvars,tsoftw,TMP,decrease=False)
          dcoefs,pathvars       = calc_coefs(target,tcommon,drst,pathvars,ltemp,plotfile)
          del drst
    # print summary with the coefficients
    fncs.print_string(PS.spath_allcoefs(ltemp,dcoefs),3)
    # return data
    return dcoefs, pathvars
#===============================================================#

#===============================================================#
def execute_pfn(itargets,dchem,idata,status,case):
    reactions = set([])
    for target in itargets:
        ctc1, itc1 = PN.name2data(target)
        the_reaction = None
        for reaction,(Rs,TS,Ps) in dchem.items():
            try   : ctc2, itc2 = PN.name2data(TS)
            except: continue
            if ctc1 == ctc2:
               the_reaction = reaction
               if itc1 == itc2: break
        if the_reaction is not None: reactions.add(the_reaction)
    if len(reactions) == 0: return
    reactions = sorted(list(reactions))

    print "   The selected transitions states are involved in defined reactions!"
    pfn_targets = set([])
    for reaction in reactions:
        print "      * %s"%reaction
        Rs,TS,Ps = dchem[reaction]
        for xx in Rs+[TS]+Ps: pfn_targets.add(PN.name2data(xx)[0])
    if len(pfn_targets) == 0: return
    pfn_targets = sorted(list(pfn_targets))
    print

    import opt_pfn as pfn
    print "   Calculating partition functions for the next targets:"
    for target in pfn_targets:
        print "      * %s"%target
    print
    pfn.main(idata,status,case,targets=pfn_targets)
#===============================================================#

#===============================================================#
def main(idata,status,case,targets="*"):

    stat2check = [1,2,3,4]
    mustexist  = [PN.DIR1]
    tocreate   = [PN.DIR4,PN.DIR2,PN.DIR3,PN.DIR5,PN.DIR6,PN.TMP]

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


    # print selected software
    print "   Selected software: %s"%software
    print

    # read high level file
    if dlevel: dhighlvl = RW.read_highlevelfile(hlf)
    else     : dhighlvl = {}

    # update targets
    if targets == [] or targets == "*": targets = dpath.keys()
    itargets = get_itargets(targets,dpath,dctc)

    # Print targets
    if len(itargets) != 0:
       print "   The paths will be calculated for:"
       for idx in range(0,len(itargets),4):
           print "       %s"%(", ".join(itargets[idx:idx+4]))
       print
    else:
       print "   No target fits for the selected options..."
       print

    # check if TSs are involved in any reaction
    try   : execute_pfn(itargets,dchem,idata,status,case)
    except: pass

    # loop over each target
   #import os
   #import psutil
   #process = psutil.Process(os.getpid())
    for target in itargets:
        #print '==>', process.memory_info().rss/1e6
        # collect garbage
        gc.collect()
        # Get ctc and itc from name
        ctc, itc = PN.name2data(target)
        # initialize pathvars
        pathvars = copy.deepcopy(dpath[ctc])
        # isotopic modification?
        masses = get_masses(target,dctc,dimasses)
        pathvars.set_masses(masses)
        # calculate path + coefficients
        idata = (target,dlevel,software,ltemp,dctc,pathvars,dtesLL,dchem,dhighlvl,dimasses)
        try:
           dcoefs, pathvars = deal_with_path(*idata)
        except Exception as exception:
           deal_with_exception(exception)
           sys.stdout = Logger(None)
           continue
        # update file (read, check ltemp, update, write)
        sys.stdout = Logger(None)
        print "   Updating data file: %s"%dof
        dall = RW.read_alldata(dof,ltemp)[0]
        for coef in dcoefs.keys():
            dall[coef][target] = dcoefs[coef]
        RW.write_alldata(dof,ltemp,dall)
        print
        # collect garbage
        #print '==>', process.memory_info().rss/1e6
        gc.collect()
    # print WARNINGS
    if len(WARNINGS) != 0:
       print "   WARNINGS:"
       for warning in WARNINGS:
          print "      *",warning
       print
#===============================================================#
