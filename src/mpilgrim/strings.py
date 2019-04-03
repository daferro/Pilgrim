#!/usr/bin/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.strings    |
| Last Update:  2019/01/23 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*

'''

#---------------------------------------------------------------#
import os
import time
import numpy             as np
#---------------------------------------------------------------#
import names             as PN
#---------------------------------------------------------------#
import common.fncs       as fncs
import common.steepdesc  as sd
import common.internal   as intl
import common.partfns    as partfns
#---------------------------------------------------------------#
from   common.physcons   import AMU
from   common.physcons   import KB
from   common.physcons   import KCALMOL
from   common.physcons   import ML
from   common.physcons   import SECOND
from   common.criteria   import EPS_AMU
#---------------------------------------------------------------#

rckey2str = {}
rckey2str["tst"   ] = "TST"
rckey2str["cvt"   ] = "CVT"
rckey2str["tstzct"] = "TST/ZCT"
rckey2str["tstsct"] = "TST/SCT"
rckey2str["cvtzct"] = "CVT/ZCT"
rckey2str["cvtsct"] = "CVT/SCT"

sortedrcs = ["tst","cvt","tstzct","cvtzct","tstsct","cvtsct"]

#===============================================================#
#                STRINGS USED IN --pfn                          #
#===============================================================#
def title_pfn(ctc,pof):
    date  = time.strftime("%Y/%m/%d %H:%M")
    title = "Analysis of STRUC: %s"%ctc
    # the string
    string  = ""
    string += "-"*(len(title)+2)+"\n"
    string += " %s \n"%title
    string += "-"*(len(title)+2)+"\n"
    string += "\n"
    string += "    Current date: %s\n"%date
    string += "\n"
    string += "    Pilgrim output file: %s\n"%pof
    return string
#---------------------------------------------------------------#
def getstring_pfn1(ltemp,pf_tot,pf_PIB,pf_RR,pf_HO,pf_ele):
    cols = ["T (K)","tr","rot","vib","el","tot"]
    string  = ""
    string += "-".join(["-"*12     for col in cols])+"\n"
    string += "|".join([fncs.fill_string(col,12) for col in cols])+"\n"
    string += "-".join(["-"*12     for col in cols])+"\n"
    for idx,T in enumerate(ltemp):
        row = pf_PIB[idx], pf_RR[idx], pf_HO[idx], pf_ele[idx], pf_tot[idx]
        row = [" %10.2f "%T]+[" %10.3E "%value for value in row]
        string += "|".join(row)+"\n"
    string += "-".join(["-"*12     for col in cols])+"\n"
    string += "  tr : translational partition funcion per volume unit (in au)\n"
    string += "  rot: rotational partition funcion (rigid-rotor)\n"
    string += "  vib: vibrational partition funcion (harmonic-oscillator)\n"
    string += "  el : electronic partition funcion\n"
    string += "  tot: total partition funcion per unit volume (in au)\n"
    string += "  \n"
    string += "  both rot and tot include rotational symmetry number\n"
    return string
#---------------------------------------------------------------#
def getstring_pfn2(ltemp,pf_tot):
    cols = ["T (K)","tot"]
    string  = ""
    string += "-".join(["-"*12     for col in cols])+"\n"
    string += "|".join([fncs.fill_string(col,12) for col in cols])+"\n"
    string += "-".join(["-"*12     for col in cols])+"\n"
    for idx,T in enumerate(ltemp):
        row = [ " %10.2f "%T," %10.3E "%pf_tot[idx] ]
        string += "|".join(row)+"\n"
    string += "-".join(["-"*12     for col in cols])+"\n"
    return string
#---------------------------------------------------------------#
def getstring_ctc(sptype,molecules,itcs,V0,V1):
    nconfs = sum([weight for name,weight in itcs])
    string  = "Number of conformers: %i\n"%nconfs
    string += "\n"
    string += "   min(V0) = %.8f hartree (V0 = electronic energy)\n"%V0
    string += "   min(V1) = %.8f hartree (V1 = electronic energy + ZPE)\n"%V1
    string += "\n"
    string += "   Relative energies (in kcal/mol):\n"
    # table of ctc
    lls      = [4,10,10,7,10,6,4]
    cols     = ["name","V0-min(V0)","V1-min(V1)","ZPE","mass (amu)","weight","PG"]
    if sptype == 1:
       cols.append("imag.freq.")
       lls.append(10)
    head     = " | ".join( [fncs.fill_string(col,ll) for ll,col in zip(lls,cols)] )
    division = "-"*len(head)
    # begin table
    string += "   "+division+"\n"
    string += "   "+head+"\n"
    string += "   "+division+"\n"
    refmass = None
    bool_wrongmass = False
    for molecule,itc in zip(molecules,itcs):
        name,weight = itc
        weight = "%i"%weight
        V0i  = "%7.2f"%((molecule._V0  -V0)*KCALMOL)
        V1i  = "%7.2f"%((molecule._ccV1-V1)*KCALMOL)
        zpe  = "%7.2f"%((molecule._ccV1 - molecule._V0)*KCALMOL)
        mass = "%7.2f"%( molecule._mass*AMU)
        pg   = molecule._pgroup
        row  = [name, V0i, V1i, zpe, mass, weight, pg]
        if sptype == 1:
           row.append( "%7.2f"%fncs.afreq2cm(molecule._ccfreqs[0]) )
        line = " | ".join( [fncs.fill_string(val,ll) for ll,val in zip(lls,row)] )
        string += "   "+line+"\n"
        if refmass is None: refmass = molecule._mass
        diff = abs(refmass-molecule._mass)*AMU
        if diff > EPS_AMU: bool_wrongmass = True
    string += "   "+division+"\n"
    if bool_wrongmass:
       string += "\n"
       string += "   ERROR: Mass differs from one conformer to another!\n"
    return string
#---------------------------------------------------------------#
def string_contributions(itcs,ltemp,dchi):
    # get contributions
    nn = 8
    string = ""
    for idx in range(0,len(itcs),nn):
        targets = itcs[idx:idx+nn]
        line = [" T (K) "]+["%5s"%itc for itc,weight in targets]
        head = " "+" | ".join(line)+" "
        divi = "-"*len(head)
        string += divi+"\n"
        string += head+"\n"
        string += divi+"\n"
        for idx,T in enumerate(ltemp):
            line    =  ["%7.2f"%T]+["%5.3f"%dchi[itc][idx] for itc,weight in targets]
            string += " "+" | ".join(line)+" "+"\n"
        string += divi+"\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --path: for MEP                #
#===============================================================#
def smep_title(target,pathvars,pof):
    date  = time.strftime("%Y/%m/%d %H:%M")

    # title
    title = "Minimum Energy Path for %s"%target
    divis = "="*len(title)
    # begin string
    string  = divis+"\n"
    string += title+"\n"
    string += divis+"\n"
    string += "\n"
    string += "  Current date: %s\n"%date
    string += "\n"
    string += "  Pilgrim output file: %s\n"%pof
    string += "\n"
    # reference energies
    if pathvars._reaction is not None:
       string += "  This transition state IS INVOLVED in reaction: '%s'\n"%pathvars._reaction
    else:
       string += "  This transition state IS NOT INVOLVED in any reaction\n"
    if pathvars._beyondmep:
       ctc, itc = PN.name2data(target)
       # E ref
       Eref = pathvars._eref
       string += "     Eref: %.6f hartree\n"%Eref
      ## E0
      #E0   = pathvars._e0
      #if E0 is not None:
      #   E0rel= (E0-Eref)*KCALMOL
      #   string += "     E0  : %.6f hartree ==> %.4f kcal/mol\n"%(E0,E0rel)
    else:
       string += "  WARNING: unable to get Eref from any reaction...\n"
    return string
#---------------------------------------------------------------#
def smep_init(target,software,PARALLEL,var_first,var_sdbw,var_sdfw):
    # split name into info
    ctc, itc = PN.name2data(target)
    # variables first step
    ds,mu,d3,idir = var_first
    string  = "Variables for first step\n"
    string += "   ds        %.4f bohr\n"%ds
    string += "   mu        %.4f amu\n"%(mu*AMU)
    if d3 in [None,False,"no"]: string += "   cubic     no\n"
    else                      : string += "   cubic     %s \n"%d3
    string += "   idir      %s %s\n"%idir
    string += "\n"
    # variables MEP
    method,mu,ds,sbw,hsteps,epse,epsg = var_sdbw
    method,mu,ds,sfw,hsteps,epse,epsg = var_sdfw
    string += "Variables for steepest descent\n"
    string += "   method    %s\n"%method
    string += "   mu        %.4f amu\n"%(mu*AMU)
    string += "   ds        %.4f bohr\n"%ds
    string += "   sbw       %+.4f bohr\n"%sbw
    string += "   sfw       %+.4f bohr\n"%sfw
    string += "   hsteps    %i\n"%hsteps
    string += "   epse      %.2E hartree\n"%epse
    string += "   epsg      %.2E hartree/bohr\n"%epsg
    string += "\n"
    # software and parallization
    string += "Other variables:\n"
    string += "   software  %s (for single-point calcs)\n"%software
    if PARALLEL:
       string += "   paral     yes (both sides of MEP will be calculated in unison)\n"
    else:
       string += "   paral     no  (both sides of MEP will NOT be calculated in unison)\n"
    return string
#---------------------------------------------------------------#
def smep_ff(f1,f2,f3,rstfile,xyzfile):
    string  = "Folders of interest:\n"
    string += "   temporal for calcs: %s\n"%f1
    string += "   folder with rst files : %s\n"%f2
    string += "   folder with xyz files : %s\n"%f3
    string += "\n"
    string += "Files of interest:\n"
    string += "   rst file: %s\n"%rstfile
    string += "   xyz file: %s\n"%xyzfile
    return string
#---------------------------------------------------------------#
def smep_rst(rstfile,drst):
    lbw, lfw, sbw, sfw, Ebw, Efw = sd.rstlimits(drst)
    if not os.path.exists(rstfile): return ""
    string  = "Restart file (rst) found!\n"
    if sbw is None:
        string += "   no data in it :(\n"
    else:
        ml = max(len(lbw),len(lfw))
        lbw = "%%-%is"%ml%lbw
        lfw = "%%-%is"%ml%lfw
        string += "   first point (%s): %+7.3f bohr (%+11.6f hartree)\n"%(lbw,sbw,Ebw)
        string += "   last  point (%s): %+7.3f bohr (%+11.6f hartree)\n"%(lfw,sfw,Efw)
    return string
#---------------------------------------------------------------#
def smep_ts(ts):
    ''' ts is an instance of Molecule'''
    string  = "\n"
    string += "Information about transition state:\n"
    string += "\n"
    string += "   File with info    : %s\n"%ts._gts
    string += ts.info_string(3)
    return string
#---------------------------------------------------------------#
def smep_first(symbols,xms,v0,v1):
    natoms  = fncs.howmanyatoms(v0)
    string  = "coordinates in mass-scaled (in bohr):\n"
    for at,symbol in enumerate(symbols):
        vx,vy,vz = fncs.xyz(xms,at)
        string += "  %2s   %+10.6f  %+10.6f  %+10.6f\n"%(symbol,vx,vy,vz)
    string += "\n"
    string += "v0 vector (mass-scaled):\n"
    for at,symbol in enumerate(symbols):
        vx,vy,vz = fncs.xyz(v0,at)
        string += "  %2s   %+10.6f  %+10.6f  %+10.6f\n"%(symbol,vx,vy,vz)
    if v1 is not None:
       string += "\n"
       string += "v1 vector (mass-scaled):\n"
       for at,symbol in enumerate(symbols):
           vx,vy,vz = fncs.xyz(v1,at)
           string += "  %2s   %+10.6f  %+10.6f  %+10.6f\n"%(symbol,vx,vy,vz)
    return string
#---------------------------------------------------------------#
def smep_table(drst,Eref):
    string  = "Reference energy (Eref) set to: %.6f hartree\n"%Eref
    string += "\n"
    string += "     -------------------------------------------------------------------\n"
    string += "      s (bohr) |   E (hartree)   | E-Eref (hartree) | E-Eref (kcal/mol) \n"
    string += "     -------------------------------------------------------------------\n"
    points = sd.sorted_points(drst,hess=True)
    for point in points:
        s_i, E_i = drst[point][0:2]
        dE_i      = (E_i-Eref)
        dE_i_kcal = (E_i-Eref)*KCALMOL
        string += "      %+8.4f | %+15.7f | %+16.7f |    %+12.3f   \n"%(s_i,E_i,dE_i,dE_i_kcal)
    string += "     -------------------------------------------------------------------\n"
    return string
#---------------------------------------------------------------#
def smep_tableDLEVEL(drst,tdleveldata,Eref):
    points, xx, yyll, yyhl = tdleveldata

    string  = "Reference energy (Eref) set to: %.6f hartree\n"%Eref
    string += "\n"
    string += "     --------------------------------------------------------------------------------------\n"
    string += "      s (bohr) | E_LL  (hartree) || E_HL  (hartree) | E-Eref (hartree) | E-Eref (kcal/mol) \n"
    string += "     --------------------------------------------------------------------------------------\n"
    for idx,point in enumerate(points):
        if drst[point][4] is None: continue
        s_i, Ell_i, Ehl_i = xx[idx],yyll[idx],yyhl[idx]
        dE_i      = (Ehl_i-Eref)
        dE_i_kcal = (Ehl_i-Eref)*KCALMOL
        string += "      %+8.4f | %+15.7f || %+15.7f | %+16.7f |    %+12.3f   \n"%(s_i,Ell_i,Ehl_i,dE_i,dE_i_kcal)
    string += "     --------------------------------------------------------------------------------------\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --path: for adiab pot          #
#===============================================================#
def sadipot_ics(ics,inter):
    nepl = 4
    string = ""
    if inter == "no": return string
    string += "Internal coordinates will be used!\n"
    ics_st,ics_ab,ics_lb,ics_it,ics_pt = intl.unmerge_ics(ics)
    for tics,lics in [("st",ics_st),("ab",ics_ab),("lb",ics_lb),("it",ics_it),("pt",ics_pt)]:
        for idx in range(0,len(lics),nepl):
            string += " "*4+" ".join(["%12s"%intl.ic2string((tics,ic)) for ic in lics[idx:idx+nepl]])+"\n"
    return string
#---------------------------------------------------------------#
def sadipot_table(ls, lV1, sAG, VAG, lV0, lcc_zpe, lic_zpe, Eref):
    string  =  "-----------------------------------------------------------------------\n"
    string +=  " s (bohr) |  V_MEP  | ZPE(cc) | ZPE(ic) |   V_adi   || Eref+V_adi (au) \n"
    string +=  "-----------------------------------------------------------------------\n"
    for idx,s_i in enumerate(ls):
        s_i   = "%+7.3f"%s_i
        V0    = "%+7.3f"%(lV0[idx]*KCALMOL)
        V1    = "%+9.3f"%(lV1[idx]*KCALMOL)
        V1_au = "%.8f"%(lV1[idx]+Eref)
        cczpe = "%7.3f"%(lcc_zpe[idx] *KCALMOL)
        arrow = "   "
        try:
          iczpe = "%7.3f"%(lic_zpe[idx]*KCALMOL)
          if float(cczpe)-float(iczpe) > 0.1: arrow = "<--"
        except:
          iczpe = "   -   "
        string +=  " %7s  | %7s | %7s | %7s | %9s || %15s  %s\n"%(s_i,V0,cczpe,iczpe,V1,V1_au,arrow)
    string +=  "-----------------------------------------------------------------------\n"
    # Maximum
    sAG = "%+7.3f"%sAG
    VAG_au = "%.8f"%(VAG+Eref)
    VAG = "%+9.3f"%(VAG*KCALMOL)
    string +=  " %7s  |      (maximum in Vadi)      | %9s || %15s\n"%(sAG,VAG,VAG_au)
    string +=  "-----------------------------------------------------------------------\n"
    string += "\n"
    return string
#---------------------------------------------------------------#
def sadipot_freqs(ls,lcc_frqs,lic_frqs):
    nepl = 5
    # Generate matrix with all data
    matrix_cc = fncs.ll2matrix(lcc_frqs,varcomp=None,pos=0)
    matrix_cc = [list(xxx) for xxx in zip(*matrix_cc)]
    nrows, ncols = len(matrix_cc), len(matrix_cc[0])
    # same for internal
    matrix_ic = fncs.ll2matrix(lic_frqs,varcomp=None,pos=0)
    if matrix_ic is not None: matrix_ic = [list(xxx) for xxx in zip(*matrix_ic)]
    # Print frequencies
    string = ""
    for idxi in range(0,ncols,nepl):
        svalues = " idx | "+"|".join(["   %+8.4f   "%ss for ss in ls[idxi:idxi+nepl]])
        string +=  "-"*len(svalues)+"\n"
        string +=          svalues +"\n"
        string +=  "-"*len(svalues)+"\n"
        # add all frequencies
        for row in range(nrows):
            string += " %3i | "%row
            for col in range(idxi,min(idxi+nepl,ncols)):
                fcc = matrix_cc[row][col]
                if fcc is not None: string +=  "%5.0f "%fncs.afreq2cm(fcc)
                else              : string +=  "  -   "
                # internal
                if matrix_ic is None: fic = None
                else                : fic = matrix_ic[row][col]
                if fic is not None: string += "(%5.0f)"%fncs.afreq2cm(fic)
                else              : string += "(  -  )"
                if col+1 == min(idxi+nepl,ncols): string += "\n"
                else                            : string += " |"
    return string
#---------------------------------------------------------------#
def sadipot_checks(checks):
    ok1, ok2, ok3 = checks
    string = "Comparing Cartesian coordinates (cc) vs internal coordinates (ic) frequencies...\n"
    if not ok1: string += "   * There are points with different numbers of frequencies!\n"
    if not ok2: string += "   * Data at transition state differ!\n"
    if not ok3: string += "   * ZPE(cc) - ZPE(ic) > 0.1 kcal/mol at some point(s)!\n"
    if not ok1 or not ok2:
       string += "   * cc adiabatic potential was chosen!\n"
    if False     in checks:
        string += "   * WARNING! MAYBE THE RESULTS ARE NOT RELIABLE! :/\n"
    if False not in checks:
        string += "   * Everything seems fine! :)\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --path: for CVT                #
#===============================================================#
def scvt_gibbs(svals,temps,mgibbs,pathvars,gibbsTS):
    #---------------------------------------
    def get_table_cvtgibbs(svals,temps,mgibbs):
        nrows = len(svals)
        nepl = 5
        matrix_string = ""
        for idx in range(0,len(temps),nepl):
            colnames = ["  s value  "]+[" %7.2f K "%T for T in temps[idx:idx+nepl]]
            thead    = " "+"|".join(colnames)+" "
            tdiv     = "-"*len(thead)
            matrix_string += tdiv +"\n"
            matrix_string += thead+"\n"
            matrix_string += tdiv +"\n"
            for row in range(nrows):
                values  = [" %+9.4f "%svals[row]]
                values += [" %9.3f "%(gibbs*KCALMOL) for gibbs in mgibbs[row][idx:idx+nepl]]
                line    = " "+"|".join(values)+" "
                matrix_string += line+"\n"
            matrix_string += tdiv +"\n\n"
        return matrix_string
    #---------------------------------------
    nrows, ncols = mgibbs.shape
    string = ""
    if gibbsTS is not None:
       string += "Matrix of gibbs free energy with regards to TS [kcal/mol]\n"
       string += "    - shape: %i x %i\n"%(nrows,ncols)
       string += "    - sigma_rot was set to 1 along MEP for the\n"
       string += "      calculation of Gibbs free energy profile\n"
       string += "\n"
       string += get_table_cvtgibbs(svals,temps,mgibbs)
    # gibbs matrix with regards reactants
    if pathvars._GibbsR is not None:
       DGRTS = [gTS-gR for gR,gTS in zip(pathvars._GibbsR,gibbsTS)]
       for col in range(ncols):
           for row in range(nrows):
               mgibbs[row][col] += DGRTS[col]
       string += "Matrix of gibbs free energy with regards to '%s' reactants [kcal/mol]\n"%pathvars._reaction
       string += "    - shape: %i x %i\n"%(nrows,ncols)
       string += "    - sigma_rot was set to 1 along MEP for the\n"
       string += "      calculation of Gibbs free energy profile\n"
       string += "\n"
       string += get_table_cvtgibbs(svals,temps,mgibbs)
    return string
#---------------------------------------------------------------#
def scvt_coefs(lcvt_s, lcvt_gamma, temps):
    string  = "CVT variational coefficient: \n"
    string += "    - variation in Gibbs free energy (DGFE) in kcal/mol\n"
    string += "            DGFE = -R T ln(Gamma_CVT) \n"
    string += "\n"
    string += "    -----------------------------------------------\n"
    string += "      T (K)  |  s_CVT  |  Gamma_CVT  |    DGFE     \n"
    string += "    -----------------------------------------------\n"
    for idx,T in enumerate(temps):
        s_i, gamma_i = lcvt_s[idx], lcvt_gamma[idx]
        gibbs        = - KB * T * np.log(gamma_i) * KCALMOL
        string += "     %7.2f | %+7.4f | %11.4E | %9.4f \n"%(T, s_i,gamma_i,gibbs)
    string += "    -----------------------------------------------\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --path: for SCT                #
#===============================================================#
def ssct_init(E0,VadiSpl,pathvars,v1mode):
    sbw, Vbw = VadiSpl.get_alpha()
    sts, Vts = VadiSpl.get_saddle()
    sfw, Vfw = VadiSpl.get_omega()
    sAG, VAG = VadiSpl.get_max()
    V1R      = pathvars._V1R
    V1P      = pathvars._V1P
    string  = "Energy limits for tunneling [E0,VAG]:\n"
    string += "\n"
    string += "  ---------------------------------------\n"
    string += "            | s (bohr) | Vadi (kcal/mol) \n"
    string += "  ---------------------------------------\n"
    if V1R is None:
       string += "  reactants |          |    --    \n"
    else:
       V1R = V1R - pathvars._eref
       string += "  reactants |          | %+8.3f \n"%(V1R*KCALMOL)
    string += "  bw        | %+8.4f | %+8.3f \n"%(sbw,Vbw*KCALMOL)
    string += "  saddle    | %+8.4f | %+8.3f \n"%(sts,Vts*KCALMOL)
    string += "  fw        | %+8.4f | %+8.3f \n"%(sfw,Vfw*KCALMOL)
    if V1P is None:
       string += "  products  |          |    --    \n"
    else:
       V1P = V1P - pathvars._eref
       string += "  products  |          | %+8.3f \n"%(V1P*KCALMOL)
    string += "  =======================================\n"
    string += "  E0        |          | %+8.3f \n"%(E0*KCALMOL)
    string += "  VAG       | %+8.4f | %+8.3f \n"%(sAG,VAG*KCALMOL)
    string += "  =======================================\n"
    string += "\n"
    if type(pathvars._e0) == float:
       string += "  E0 defined by user in '%s'\n"%PN.IFILE3
       string += "\n"
    string += "v1mode  : %s\n"%v1mode
    string += "muintrpl: %s %i\n"%(pathvars._muintrpl)
    return string
#---------------------------------------------------------------#
def ssct_mueff(svals,VadiSpl,lkappa,ltbar,lmueff,toignore=[]):
    string  = "Summary table\n"
    string += "  - s in bohr\n"
    string += "  - kappa (curvature) in bohr^-1\n"
    string += "  - Turning point in bohr\n"
    string += "  - Effective mass (mueff) in a.u.\n"
    string += "  ---------------------------------------------------------\n"
    string += "      s    |   V_adi   |   kappa   | turnpoint |  mueff    \n"
    string += "  ---------------------------------------------------------\n"
    for idx in range(len(svals)):
         s      = "%+7.3f"%svals[idx]
         Vadi   = "%+9.3f"%(VadiSpl(svals[idx])*KCALMOL)
         kappa  = "%+9.2E"%lkappa[idx]
         tbar   = "%9.5f"%ltbar[idx]
         mu_eff = "%9.4f"%(lmueff[idx])
         string += "   %s | %s | %s | %s | %s "%(s,Vadi,kappa,tbar,mu_eff)
         if idx in toignore: string += "**"
         string += "\n"
    string += "  ---------------------------------------------------------\n"
    string += "   ** kappa, turnpoint and mueff were interpolated\n"
    return string
#---------------------------------------------------------------#
def ssct_E0VAG(E0,VAG):
    string  = "Transmission probabilities will be calculated between E0 and VAG:\n"
    string += "   E_0    = %8.4f kcal/mol\n"%( E0 *KCALMOL )
    string += "   V^{AG} = %8.4f kcal/mol\n"%( VAG*KCALMOL )
    string += "\n"
    return string
#---------------------------------------------------------------#
def ssct_probs(E_list,probs_ZCT,probs_SCT,rpoints_SCT):
    string = ""
    # Get P(E) for all these energies
    head = " E [kcal/mol]  |  P^ZCT(E)  |  P^SCT(E)  |  Integration limits [bohr] "
    divi = "-"*len(head)
    string += "  "+divi+"\n"
    string += "  "+head+"\n"
    string += "  "+divi+"\n"
    for idx in range(len(E_list)):
        E           = "%10.4f"%(E_list[idx]*KCALMOL)
        pE_ZCT      = "%10.3e"%(probs_ZCT[idx])
        pE_SCT      = "%10.3e"%(probs_SCT[idx])
        rpoints     = rpoints_SCT[idx]
        # string return points
        rp_string = ""
        for idx,(start,end) in enumerate(rpoints):
            rp_string += "[%+.3f,%+.3f]U"%(start,end)
            if (idx+1)%2==0 and idx+1 != len(rpoints): rp_string+= "\n"+" "*46
        string += "     %s  | %s | %s |  %s\n"%(E,pE_ZCT,pE_SCT,rp_string[:-1])
    string += "  "+divi+"\n"
   #string += "   Note: Sometimes, the value of V at a given s may be greater than E.  \n"
   #string += "         In such cases, a letter may be found next to the interval.     \n"
   #string += "           * an 'L' indicates that this happens at the left-side of s.  \n"
   #string += "           * an 'R' indicates that this happens at the right-side of s. \n"
   #string += "           * an 'B' indicates that this happens at both sides.          \n"
    return string
#---------------------------------------------------------------#
def ssct_diffs(lE_SCT,diffs_SCT):
    maxdiff = 10.0
    nepl    = 5
    bigdiff = [E for E,diff in zip(lE_SCT,diffs_SCT) if diff > maxdiff]
    string = ""
    if len(bigdiff) > 0:
       string += "Action integral (theta) is calculated by:\n"
       string += "   * Gaussian quadrature (theta_g)\n"
       string += "   * trapezoidal rule    (theta_t)\n"
       string += "and theta = min(theta_g,theta_t)\n"
       string += "\n"
       string += "WARNING: differences > %i%% for the following energies:\n"%maxdiff
       for idx in range(0,len(bigdiff),nepl):
           string += "         "
           string += " ".join(["%8.4f"%(E*KCALMOL) for E in bigdiff[idx:idx+nepl]])
           string += "\n"
      #string += "   * theta = theta_t in those cases\n"
       string += "\n"
    return string

#---------------------------------------------------------------#
def ssct_kappa(temps,KAPPA,lIi,RTE,E0,case="sct"):
    string  = "%s transmission coefficient (KAPPA):\n"%(case.upper())
    string += "--------------------------------------------------------------\n"
    string += "   T (K)   |  %%I1   |  %%I2   |  %%I3   |  KAPPA_%3s  |   RTE   \n"%(case.upper())
    string += "--------------------------------------------------------------\n"
    for idx,T in enumerate(temps):
        I1, I2, I3 = lIi[idx]
        kappa, rte = KAPPA[idx], RTE[idx]*KCALMOL
        I1 *= 100./kappa
        I2 *= 100./kappa
        I3 *= 100./kappa
        string += " %9.2f | %6.2f | %6.2f | %6.2f | %11.4E | %7.3f"%(T,I1,I2,I3,kappa,rte)
        if (rte-E0*KCALMOL) < 1.5: string += " <-- RTE close to E0"
        string += "\n"
    string += "--------------------------------------------------------------\n"
    string += " RTE stands for Representative Tunneling Energy               \n"
    string += " %I1: contribution to KAPPA from energies in [E0     ,VAG    )\n"
    string += " %I2: contribution to KAPPA from energies in [VAG    ,2VAG-E0)\n"
    string += " %I3: contribution to KAPPA from energies in [2VAG-E0,infty  )\n"
    return string
#---------------------------------------------------------------#
def ssct_convergence(convlist_sct,convlist_lim,scterr):
    string  = "Convergence of kappa(SCT) for the lowest temperature \n"
    string += "-------------------------------------------------\n"
    string += "\n"
    string += "     scterr = %.2f %%\n"%scterr
    string += "\n"
    string += "     --------------------------------------------------------------------\n"
    string += "       step   |   sbw   |   sfw   | KAPPA (SCT) |  diff(%)  | converged? \n"
    string += "     --------------------------------------------------------------------\n"
    for idx, SCT in enumerate(convlist_sct):
        sbw, sfw = convlist_lim[idx]
        if idx == 0:
            diff  = "    -    "
            sconv = "    -    "
        else:
            SCT_pr = convlist_sct[idx-1]
            diff = "%9.2f"%float(100*abs(SCT-SCT_pr)/SCT)
            if float(diff) > scterr: sconv = "    NO    "
            else                   : sconv = "    YES   "
        string += "          %3i | %+7.3f | %+7.3f | %+11.4E | %s | %s \n"%(idx-1,sbw,sfw,SCT,diff,sconv)
    string += "     --------------------------------------------------------------------\n"
    # add square around for visualization
    with_square  = "="*82+"\n"
    for line in string.split("\n"):
        with_square += "|| " + "%-77s"%line + "||\n"
    with_square += "="*82+"\n"
    string = with_square
    # return string
    return string
#===============================================================#

#===============================================================#
#                STRINGS USED IN --path: for CAG                #
#===============================================================#
def scag_table(ltemp, dE_cagtst, cagtst, dE_cagcvt, cagcvt):
    string  = "    ----------------------------------------------------------------\n"
    string += "      T (K)  |  VAG-V(TS)  |  CAG_TST  || VAG-V(sCVT) |  CAG_CVT  \n"
    string += "    ----------------------------------------------------------------\n"
    for idx,T in enumerate(ltemp):
        dE1 = "%11.4f"%(dE_cagtst[idx]*KCALMOL)
        cf1 = "%9.3E"%(cagtst[idx])
        if dE_cagcvt is None or cagcvt is None:
           dE2 = "     -     "
           cf2 = "    -    "
        else:
           dE2 = "%11.4f"%(dE_cagcvt[idx]*KCALMOL)
           cf2 = "%9.3E"%(cagcvt[idx])
        string += "     %7.2f | %s | %s || %s | %s\n"%(T,dE1,cf1,dE2,cf2)
    string += "    ----------------------------------------------------------------\n"
    string += "      Energy difference, VAG-V(s_i), in kcal/mol\n"
    return string
#===============================================================#


#===============================================================#
#                                                               #
#===============================================================#
def spath_allcoefs(ltemp,dcoefs):
    coefs = ["CVT","ZCT","SCT","CAGTST","CAGCVT"]
    # print coefs again
    table = ""
    thead = [" T  (K) "]+[fncs.fill_string(coef,13) for coef in coefs]
    thead = " "+" | ".join(thead)+" "
    tdivi = "-"*len(thead)
    table += tdivi+"\n"
    table += thead+"\n"
    table += tdivi+"\n"
    for idx,T in enumerate(ltemp):
        tline = ["%8.2f"%T]
        for coef in coefs:
            coef = coef.lower()
            if coef in dcoefs.keys():
                value = "%13.4E"%dcoefs[coef][idx]
            else:
                value = "      -      "
            tline.append(value)
        tline = " "+" | ".join(tline)+" "
        table += tline+"\n"
    table += tdivi+"\n"
    # return string
    string  = "SUMMARY OF CALCULATED COEFFICIENTS:\n"
    string += "\n"
    for line in table.split("\n"): string += "    "+line+"\n"
    string += "\n"
    return string
#===============================================================#


#===============================================================#
#                STRINGS USED IN --rcons                        #
#===============================================================#
def title_rcons(rcname,pof):
    date  = time.strftime("%Y/%m/%d %H:%M")
    # Title
    title    = "| Reaction to analyze: %s |"%rcname
    division = "-"*len(title)
    # the string
    string  = division+"\n"
    string += title   +"\n"
    string += division+"\n"
    string += "\n"
    string += "    Current date: %s\n"%date
    string += "\n"
    string += "    Pilgrim output file: %s\n"%pof
    return string
#---------------------------------------------------------------#
def conf_table(targets,dctc,dall,V0R,V1R):
    # for not defined TS
    if None in targets: return []
    # get key for each conformer
    allkeys = []
    for target in targets:
        keys = []
        # get ctc and itc
        if "." in target: ctc,itc = target.split(".")
        else            : ctc,itc = target, None
        # ctc in dctc?
        if ctc not in dctc.keys(): return []
        # get list of itcs
        if itc is None: itcs = dctc[target]._itcs
        else          : itcs = [(itc,1)]
        # gen keys
        keys = [PN.struckey(ctc,itc) for itc,weight in itcs]
        # append data
        allkeys.append(keys)
    # only if len == 1 or 2
    if   len(allkeys) == 1: allkeys = allkeys[0]
    elif len(allkeys) == 2:
       keys1, keys2 = allkeys
       allkeys = [k1+" + "+k2 for k1 in keys1 for k2 in keys2]
    else: return []
    # generate lines of table
    lines = []
    for thekey in allkeys:
        keys = thekey.split("+")
        V0i,V1i = 0.0, 0.0
        for ki in keys:
            ctc,itc = ki.split(".")
            ctc     = ctc.strip()
            itc     = itc.strip()
            try   : V0,V1 = dall["pfn"][PN.struckey(ctc,itc)][0:2]
            except: return []
            V0i += V0
            V1i += V1
        line = (thekey, (V0i-V0R), (V1i-V1R))
       #line = (thekey, (V0i-V0R), (V1i-V1R), (V1i-V0R))
        lines.append(line)
    return lines
#---------------------------------------------------------------#
def srcons_anhtable(trganh,ltemp,ANH4KEQ,ANH4RCF,ANH4RCB,bools):
    iblank = "    "
    string  = "Anharmonicity ratios were found for:\n"
    string += "\n"
    for trg in trganh: string += iblank+"* %s\n"%trg
    string += "\n"
    string += "Anharmonic corrections (ANHC) for equilibrium (Keq) and rate (k) constants:\n"
    string += "\n"
    thead = "  T  (K)  |  ANHC(Keq)  |  ANHC(k,fw)  |  ANHC(k,bw)  "
    tdivi = "-"*len(thead)
    string += iblank+tdivi+"\n"
    string += iblank+thead+"\n"
    string += iblank+tdivi+"\n"
    for idx,T in enumerate(ltemp):
        if bools[0]: col1 = "%11.3E"%ANH4KEQ[idx]
        else       : col1 = "     -     "
        if bools[1]: col2 = "%11.3E"%ANH4RCF[idx]
        else       : col2 = "     -     "
        if bools[2]: col3 = "%11.3E"%ANH4RCB[idx]
        else       : col3 = "     -     "
        tline = " %8.2f | %s | %s | %s "%(T,col1,col2,col3)
        string += iblank+tline+"\n"
    string += iblank+tdivi+"\n"
    return string
#---------------------------------------------------------------#
def srcons_relenergies(dctc,dall,Rs,TS,Ps,V0R,V1R):
    dall["pfn"]
    lines = []
    lines += conf_table( Rs ,dctc,dall,V0R,V1R)+[None]
    lines += conf_table([TS],dctc,dall,V0R,V1R)+[None]
    lines += conf_table( Ps ,dctc,dall,V0R,V1R)+[None]
    ml   = max([len(line[0]) for line in lines if line is not None])
    col1 =("%%-%is"%ml)%"i"
    thead = " %s | V0(i)-V0 | V1(i)-V1 "%col1
   #thead = " %s | V0(i)-V0 | V1(i)-V1 | V1(i)-V0 "%col1
    tdivi = "-"*len(thead)
    string  = "Relative energies (kcal/mol):\n"
    string += "\n"
    string += "   "+"V0(i) is the electronic energy of 'i'\n"
    string += "   "+"V1(i) is the electronic energy + the harmonic oscillator zero point energy (ZPE) of 'i'\n"
    string += "   "+"V1(i) = V0(i)+ZPE(i)\n"
    string += "\n"
    string += "   "+"min{V0(i)} of reactants ==> V0 = %.8f hartree\n"%V0R
    string += "   "+"min{V1(i)} of reactants ==> V1 = %.8f hartree\n"%V1R
    string += "\n"
    string += "   "+tdivi+"\n"
    string += "   "+thead+"\n"
    string += "   "+tdivi+"\n"
    for idx,line in enumerate(lines):
        if line is None:
           if idx+1 == len(lines): continue
           if lines[idx-1] is None and lines[idx+1] is None: continue
           if lines[idx+1] is None: continue
           string += "   "+tdivi+"\n"
           continue
        target, V00i, V11i = line
       #target, V00i, V11i, V01i = line
        V00i *= KCALMOL
        V11i *= KCALMOL
        #V01i *= KCALMOL
        target = ("%%-%is"%ml)%target
        tline  = " %s | %8.2f | %8.2f "%(target,V00i, V11i)
       #tline  = " %s | %8.2f | %8.2f | %8.2f"%(target,V00i, V11i, V01i)
        string += "   "+tline+"\n"
    string += "   "+tdivi+"\n"
    string += "\n"
    return string
#---------------------------------------------------------------#
def srcons_keq(ltemp,keqs,nR,nP,wfw=1.0,wbw=1.0,case=0):
    # calculate Gibbs
    gibbs = partfns.Kc2GFE(ltemp,keqs,dn=nR-nP)
    # begin table
    if case == 0: string = "Harmonic Equilibrium Constants:\n"
    else        : string = "Anharmonic Equilibrium Constants:\n"
    string += "\n"
    if nR == nP: string += "  - Units of Keq (R-->P): dimensionless\n"
    else       : string += "  - Units of Keq (R-->P): (molecule/mL)^%i\n"%(nR-nP)
    string += "  - Units of Gibbs: kcal/mol\n"
    if wfw!=1.0 or wbw!=1.0:
        string += "  - Equilibrium constants are corrected by a factor of: %i/%i\n"%(wfw,wbw)
    string += "\n"
    human_units = ML**(nR-nP)
    string += "   ---------------------------------------------------\n"
    string += "     T (K)  | Gibbs (R-->P) | Keq(R-->P) | Keq(P-->R) \n"
    string += "   ---------------------------------------------------\n"
    for idx,T in enumerate(ltemp):
        keq   = keqs[idx]  * human_units
        dG    = gibbs[idx] * KCALMOL
        string += "    %7.2f |  %+11.3f  | %10.3E | %10.3E \n"%(T,dG,keq,1.0/keq)
    string += "   ---------------------------------------------------\n"
    return string
#---------------------------------------------------------------#
def srcons_rate(ltemp,rates,nR,rctype):
    human_units = ML**(nR-1.0) / SECOND
    gfes  = partfns.rate2GFE(ltemp,rates,dn=nR-1)
    rcstr = rckey2str[rctype]
    string  = "---------------------------------------\n"
    string += "   T (K)   |   %9s   |   Gibbs   \n"%fncs.fill_string(rcstr,9)
    string += "---------------------------------------\n"
    for idx,T in enumerate(ltemp):
        rate        = rates[idx] * human_units
        gfe         = gfes[idx] * KCALMOL
        string += " %9.2f |  %11.3E  | %9.3f \n"%(T,rate,gfe)
    string += "---------------------------------------\n"
    return string
#---------------------------------------------------------------#
def srcons_ratios(ltemp,ratios,litcs,rcs):
    nepl = 8
    string  = ""
    for rctype in rcs:
        string += "- For %s rate constant:\n"%(rckey2str[rctype])
        ibs = "    "
        for idx in range(0,len(litcs),nepl):
            head = "  T  (K)  |" + "|".join([" %5s "%tsi for tsi,weight in litcs[idx:idx+nepl]])
            string += ibs+"-"*len(head)+"\n"
            string += ibs+        head +"\n"
            string += ibs+"-"*len(head)+"\n"
            for idx2,T in enumerate(ltemp):
                string += ibs+" %8.2f |"%T
                string += "|".join([" %5.3f "%(ratios[rctype][tsi][idx2]*weight) \
                              for tsi,weight in litcs[idx:idx+nepl]])
                string += "\n"
            string += ibs+"-"*len(head)+"\n"
        string += "\n"
    return string
#---------------------------------------------------------------#
def srcons_tstratios(ltemp,tst_ratios,case="MSHO"):
    itcs = sorted(tst_ratios.keys())
    if len(itcs) == 1: return ""
    ncols = 6
    blank = "    "    

    if case in rckey2str.keys(): case = rckey2str[case]
    table  = "Contribution of each TS conformer to %s:\n"%case
    table +=  "\n"
    for idx1 in range(0,len(itcs),ncols):
        targets = itcs[idx1:idx1+ncols]
        head = " "+" | ".join([" T (K) "]+[" %3s "%target for target in targets])+" "
        division = "-"*len(head)
        table += blank+division+"\n"
        table += blank+head    +"\n"
        table += blank+division+"\n"
        for idx2,T in enumerate(ltemp):
            linedata = ["%7.2f"%T]
            for itc in targets:
                linedata.append( "%5.3f"%(tst_ratios[itc][idx2]) )
            line = " "+" | ".join(linedata) +" "
            table += blank+line+"\n"
        table += blank+division+"\n\n"
    return table
#---------------------------------------------------------------#
def srcons_vtstratios(ltemp,tst_ratios,corrfactors):
    # get name of correction coefficients
    rctypes = corrfactors.keys()
    sorted1 = [rctype for rctype in sortedrcs if rctype     in rctypes  ]
    sorted2 = [rctype for rctype in rctypes   if rctype not in sortedrcs]
    rctypes = sorted1 +  sorted(sorted2)
    # get string
    whole_string = ""
    for rctype in rctypes:
        try:
          cf_tot = corrfactors[rctype]["total"]
          dratios = {}
          for itc in tst_ratios.keys():
              cf_i   = corrfactors[rctype][itc]
              chi_i  = tst_ratios[itc]
              contri = cf_i*chi_i/cf_tot
              dratios[itc] = contri
          if len(dratios.keys()) < 2: continue
          whole_string += srcons_tstratios(ltemp,dratios,case=rctype)
          #print_string(PS.srcons_tstratios(ltemp,dratios,case=rctype),7)
        except: pass
    return whole_string
#---------------------------------------------------------------#
def string_corrfactors(ltemp,corrfactors):
    empty = True
    # get name of correction coefficients
    rctypes = corrfactors.keys()
    sorted1 = [rctype for rctype in sortedrcs if rctype     in rctypes  ]
    sorted2 = [rctype for rctype in rctypes   if rctype not in sortedrcs]
    rctypes = sorted1 +  sorted(sorted2)
    # table head and table division
    table_head = " | ".join( ["T(K) \ X"]+["%9s"%rckey2str[rctype] for rctype in rctypes] )
    table_head = " "+table_head+" "
    division   = len(table_head)*"-"
    # start table
    string  = "Total correction coefficients to TST rate constant:\n"
    string += "\n"
    string += "   "+"<gamma>^X = k^X / k^TST"+"\n"
    string += "\n"
    string += "   "+division  +"\n"
    string += "   "+table_head+"\n"
    string += "   "+division  +"\n"
    for idx,T in enumerate(ltemp):
        values = ["%9.2f"%T]
        for rctype in rctypes:
            try   :
              values += ["%9.3E"%(corrfactors[rctype]["total"][idx])]
              empty = False
            except: values += ["    -    "]
        string += "   "+" | ".join(values)+"\n"
    string += "   "+division+"\n"
    if empty: return ""
    else    : return string
#----------------------------------------------------------#
def string_rcons(ltemp,rcons,nR,weight=1.0,case=0):
    human_units = ML**(nR-1.0) / SECOND
    # get name of correction coefficients
    if case == 0:
       rctypes = [rctype for rctype in rcons.keys() if not rctype.startswith("anh")]
    else:
       rctypes = [rctype for rctype in rcons.keys() if     rctype.startswith("anh")]
       rctypes = [rctype.split("anh")[1] for rctype in rctypes]
    sorted1 = [rctype for rctype in sortedrcs if rctype     in rctypes  ]
    sorted2 = [rctype for rctype in rctypes   if rctype not in sortedrcs]
    rctypes = sorted1 +  sorted(sorted2)
    # table head and table division
    table_head = " | ".join( ["  T (K)  "]+["%9s"%rckey2str[rctype] for rctype in rctypes] )
    division   = len(table_head)*"-"
    # start table
    if case == 0: string  = "Harmonic Rate Constants "
    else        : string  = "Anharmonic Rate Constants "
    if nR == 1  : string += "in sec^-1\n"
    else        : string += "in cc^%i/(molecule^%i sec)\n"%(nR-1,nR-1)
    string += "\n"
    if weight != 1.0:
       string += "   - Weight for the rate constant is applied: %.3f\n"%weight
       string += "\n"
    string += division  +"\n"
    string += table_head+"\n"
    string += division  +"\n"
    for idx,T in enumerate(ltemp):
        values = ["%9.2f"%T]
        for rctype in rctypes:
            if case == 0: key = rctype
            else        : key = "anh"+rctype
            try   : values += ["%9.3E"%(rcons[key][idx]*human_units)]
            except: values += ["    -    "]
        string += " | ".join(values)+"\n"
    string += division+"\n"
    string += "\n"
    return string
#----------------------------------------------------------#
def string_gibbsactenergy(ltemp,rcons,nR,case=0):
    # get name of correction coefficients
    if case == 0:
       rctypes = [rctype for rctype in rcons.keys() if not rctype.startswith("anh")]
    else:
       rctypes = [rctype for rctype in rcons.keys() if     rctype.startswith("anh")]
       rctypes = [rctype.split("anh")[1] for rctype in rctypes]
    sorted1 = [rctype for rctype in sortedrcs if rctype     in rctypes  ]
    sorted2 = [rctype for rctype in rctypes   if rctype not in sortedrcs]
    rctypes = sorted1 +  sorted(sorted2)
    # table head and table division
    table_head = " | ".join( ["  T (K)  "]+["%9s"%rckey2str[rctype] for rctype in rctypes] )
    division   = len(table_head)*"-"
    # start table
    if case == 0: string  = "Harmonic Gibbs Free Energy of Activation in kcal/mol\n"
    else        : string  = "Anharmonic Gibbs Free Energy of Activation in kcal/mol\n"
    string += "calculated as GFEA = -kB T ln(h k / kB T) \n"
    string += "\n"
    string += division  +"\n"
    string += table_head+"\n"
    string += division  +"\n"
    for idx,T in enumerate(ltemp):
        values = ["%9.2f"%T]
        rates = rcons[rctype]
        for rctype in rctypes:
            if case == 0: key = rctype
            else        : key = "anh"+rctype
            rate = rcons[key][idx]
            gfe  = partfns.rate2GFE([T],[rate],dn=nR-1)[0]
            try   : values += ["%9.3f"%(gfe*KCALMOL)]
            except: values += ["    -    "]
        string += " | ".join(values)+"\n"
    string += division+"\n"
    string += "\n"
    return string
#===============================================================#


#===============================================================#
#                                                               #
#===============================================================#
def sfit_rcons(ltemp,rcons,nR,case=0):
    human_units = ML**(nR-1.0) / SECOND
    # get name of correction coefficients
    if case == 0:
       rctypes = [rctype for rctype in rcons.keys() if not rctype.startswith("anh")]
    else:
       rctypes = [rctype for rctype in rcons.keys() if     rctype.startswith("anh")]
       rctypes = [rctype.split("anh")[1] for rctype in rctypes]
    sorted1 = [rctype for rctype in sortedrcs if rctype     in rctypes  ]
    sorted2 = [rctype for rctype in rctypes   if rctype not in sortedrcs]
    rctypes = sorted1 +  sorted(sorted2)
    # table head and table division
    table_head = " | ".join( ["  T (K)  "]+["%9s"%rckey2str[rctype] for rctype in rctypes] )
    division   = len(table_head)*"-"
    # start table
    if case == 0: string  = "Harmonic Rate Constants "
    else        : string  = "Anharmonic Rate Constants "
    if nR == 1  : string += "in sec^-1\n"
    else        : string += "in cc^%i/(molecule^%i sec)\n"%(nR-1,nR-1)
    string += "\n"
    string += "   "+division  +"\n"
    string += "   "+table_head+"\n"
    string += "   "+division  +"\n"
    for idx,T in enumerate(ltemp):
        values = ["%9.2f"%T]
        for rctype in rctypes:
            if case == 0: key = rctype
            else        : key = "anh"+rctype
            try   : values += ["%9.3E"%(rcons[key][idx]*human_units)]
            except: values += ["    -    "]
        string += "   "+" | ".join(values)+"\n"
    string += "   "+division+"\n"
    return string
#----------------------------------------------------------#
def sfit_anafit(dfit,key):
    rctype, rcname, direc = key.split(".")
    if rctype.startswith("anh"): rctype = rctype.split("anh")[1]
    string = ""
    #-------------#
    # fitting     #
    #-------------#
    # print tables
    cols   = [1,2,3,4]
    rows   = ["A","B","n","Tr","T0"]
    rcstr  = rckey2str[rctype]
    string += "-----------------------------------------------------------------\n"
    string += " %7s |     (1)     |     (2)     |     (3)     |     (4)     \n"%rcstr
    string += "-----------------------------------------------------------------\n"
    # add parameters
    for idx,row in enumerate(rows):
        string += "   %3s   |"%row
        for col in cols:
            params,r2 = dfit.get(col,(None,None))
            try   :
                if col == 4: params = params[0:3]+[params[4],params[3]]
                param = "%11.3E"%params[idx]
            except:
                param = " "*11
            string += " %s "%param + "|"
        string = string[:-1]+"\n"
    string += "-----------------------------------------------------------------\n"
    # add 1-r^2
    string += "   1-r^2 |"
    for col in cols:
        params,r2 = dfit.get(col,(None,None))
        try   : r2 = "%7.2E"%(1-r2)
        except: r2 = " "*7
        string += " %11s |"%r2
    string = string[:-1]+"\n"
    string += "-----------------------------------------------------------------\n"
    # lines to use analytical expressions in KMC (just in case)
    lines = {}
    for anatype in dfit.keys():
        coefs, r2 = dfit[anatype]
        if anatype == 4: coefs = coefs[0:3]+[coefs[4],coefs[3]]
        coefs = " ".join(["%10.3E"%coef for coef in coefs])
        line  = "k(%s.%s) analytic%i %-55s # r^2 = %.8f\n"%(rcname,direc,anatype,coefs,r2)
        lines[anatype] = line
    return string, lines
#===============================================================#


#===============================================================#
#                STRINGS USED IN --kmc                          #
#===============================================================#
def skmc_init(ipops,POP0,excess,rcs,psteps,volume,timeunits,dof):
    # some lengths for nice format
    try   : l1 = max([len(key) for key       in rcs.keys()])
    except: l1 = 4
    try   : l2 = max([len(v1)  for (v1,v2)   in rcs.values()])
    except: l2 = 4
    try   : l3 = max([len(key) for key,val   in ipops.items() if val != 0.0])
    except: l3 = 4
    # generate string
    string  = ""
    string += "---------------------\n"
    string += "KMC basic information\n"
    string += "---------------------\n"
    string += "  Writing step set to    : %s\n"%psteps
    string += "  Simulation volume      : %.3E mL\n"%(volume*ML)
    string += "  Time units             : %s\n"%(timeunits)
    string += "\n"
    string += "  Initial populations, pop(i;t=0):\n"
    for species,ipop in sorted(ipops.items()):
        if ipop == 0.0: continue
        species = ("%%-%is"%l3)%species
        string += "      pop(%s;t=0) = %.3E particles"%(species,ipop)
        if species in excess: string += " (excess)\n"
        else                : string += "\n"
    string += "\n"
    string += "  Population of limiting reactant (POP0): %.3E particles\n"%POP0
    string += "\n"
    string += "  Rate constants:\n"
    for rname,(ktype,weight,coefs) in sorted(rcs.items()):
        if rname.endswith(".both"): rname = rname[:-5]
        rname = ("%%-%is"%l1)%rname
        ktype = ("%%-%is"%l2)%ktype
        string += "      %s for %s"%(ktype,rname)
        if weight != 1: string += " (x%i)"%weight
        string += "\n"
    string += "---------------------\n"
    string += "\n"

    return string
#---------------------------------------------------------------#
def skmc_processes(T,processes):
    string  = ""
    string += "Processes to be considered:\n"
    ml1 = max([len(" + ".join(Rs)) for Rs,Ps,k in processes])
    ml2 = max([len(" + ".join(Ps)) for Rs,Ps,k in processes])
    for Rs,Ps,k in processes:
        nR = len(Rs)
        if nR == 1:
           k *= (1.0/SECOND)
           units = "sec^-1"
        else:
            k *= (ML**(nR-1)) / SECOND
            units = "cc^%i/(molecule^%i sec)"%(nR-1,nR-1)
        str_Rs = ("%%-%is"%ml1)%" + ".join(Rs)
        str_Ps = ("%%-%is"%ml2)%" + ".join(Ps)
        string += "    %s --> %s (%.3E %s)\n"%(str_Rs,str_Ps,k,units)
    return string
#---------------------------------------------------------------#
def skmc_results(xvalues,yvalues,timeunits):
    stime    = "time (%s)"%timeunits
    nrpl     = 6
    ib       = "    "

    len1 = max([len(string) for string in yvalues.keys()+[stime]])
    molecules = sorted(yvalues.keys())

    # xvalues to units
    if timeunits == "fs" : xvalues = [t_i*SECOND*1e15  for t_i in xvalues]
    if timeunits == "ps" : xvalues = [t_i*SECOND*1e12  for t_i in xvalues]
    if timeunits == "mcs": xvalues = [t_i*SECOND*1e6   for t_i in xvalues]
    if timeunits == "ms" : xvalues = [t_i*SECOND*1e3   for t_i in xvalues]
    if timeunits == "s"  : xvalues = [t_i*SECOND       for t_i in xvalues]
    if timeunits == "min": xvalues = [t_i*SECOND/60.   for t_i in xvalues]
    if timeunits == "hr" : xvalues = [t_i*SECOND/3600. for t_i in xvalues]
    string  = ""
    string += "Evolution of population(s) with time, pop_i(t):\n"
    string += "\n"
    for idx1 in range(0,len(xvalues),nrpl):
        idx2 = idx1+nrpl
        table_head  = (" %%-%is | "%len1)%stime
        table_head += " | ".join( ["%.2E"%tt for tt in xvalues[idx1:idx2]])
        string += ib+table_head+"\n"
        string += ib+"-"*len(table_head)+"\n"
        for molecule in molecules:
            pops = yvalues[molecule][idx1:idx2]
            molecule = ("%%-%is"%len1)%molecule
            string += ib+" %s | "%molecule + " | ".join( ["%.2E"%pop for pop in pops])+"\n"
        string += "\n"
    return string
#---------------------------------------------------------------#
def skmc_finaltimes(ltemp,ftimes,timeunits):
    if timeunits == "fs" : ftimes = [t_i*SECOND*1e15  for t_i in ftimes]
    if timeunits == "ps" : ftimes = [t_i*SECOND*1e12  for t_i in ftimes]
    if timeunits == "mcs": ftimes = [t_i*SECOND*1e6   for t_i in ftimes]
    if timeunits == "ms" : ftimes = [t_i*SECOND*1e3   for t_i in ftimes]
    if timeunits == "s"  : ftimes = [t_i*SECOND       for t_i in ftimes]
    if timeunits == "min": ftimes = [t_i*SECOND/60.   for t_i in ftimes]
    if timeunits == "hr" : ftimes = [t_i*SECOND/3600. for t_i in ftimes]
    string = ""
    string += "------------------------\n"
    string += " Simulation time in %3s \n"%timeunits
    string += "------------------------\n"
    string += "\n"
    string += "     ---------------------\n"
    string += "       T (K)  | sim. time \n"
    string += "     ---------------------\n"
    for T,ftime in zip(ltemp,ftimes):
        string += "      %7.2f | %9.2E \n"%(T,ftime)
    string += "     ---------------------\n"
    string += "\n"
    return string
#---------------------------------------------------------------#
def skmc_finalratio(ltemp,species,fratios):
    NEPT    = 5 # number of elements per table
    ml      = max([len(xx) for xx in species]+[7])
    string = ""
    string += "----------------------------\n"
    string += " Final ratios (pop(i)/POP0) \n"
    string += "----------------------------\n"
    string += "\n"
    for idx in range(0,len(species),NEPT):
        head     = [" T (K) "]+[("%%%is"%ml)%xx for xx in species[idx:idx+NEPT]]
        head     = " "+" | ".join(head)+" "
        division = "-"*len(head)
        string += "     "+division+"\n"
        string += "     "+head    +"\n"
        string += "     "+division+"\n"
        for temp in ltemp:
            key1 = "%7.2f"%temp
            line = [key1]
            for species_i in species[idx:idx+NEPT]:
                value = fratios[key1][species_i]
                if value < 10: line.append( ("%%%i.3f"%ml)%value )
                else         : line.append( "%7.1E"%value )
            line = " "+" | ".join(line)+" "
            string += "     "+line+"\n"
        string += "     "+division+"\n\n"
    return string
#===============================================================#

