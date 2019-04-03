#!/usr/bin/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.pdf        |
| Last Update:  2019/01/16 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*

'''

#--------------------------------------------------#
import os
import common.steepdesc  as sd
import common.physcons   as pc
import numpy             as np
import common.files      as ff
import common.fncs       as fncs
from common.fit2anarc import anarc1
from common.fit2anarc import anarc2
from common.fit2anarc import anarc3
from common.fit2anarc import anarc4
#--------------------------------------------------#


#==================================================#
def start_end_blocks(plotfile):
    text = "".join(ff.read_file(plotfile))
    blocks = text.split("start_")[1:]
    blocks = [block.split("end_")[0] for block in blocks]
    return blocks
#--------------------------------------------------#
def cleanstr(string):
    string = string.strip()
    if string.startswith("'"): string = string[1:]
    if string.endswith("'")  : string = string[:-1]
    return string
#--------------------------------------------------#
def data_in_block(block):
    ptype   = block[0].split()[0]
    mlabel  = block[0].split()[1]
    title   = ""
    xlabel  = ""
    ylabel  = ""
    xydata  = []
    # Go line by line in block
    bool_xy = False
    for line in block[1:]:
        line = line.strip()
        if line.startswith("title" ): title  = "=".join(line.split("=")[1:])
        if line.startswith("xlabel"): xlabel = "=".join(line.split("=")[1:])
        if line.startswith("ylabel"): ylabel = "=".join(line.split("=")[1:])
        if bool_xy:
            xi, yi = line.split()
            if ptype == "plot": xx.append(float(xi))
            if ptype == "bar" : xx.append(  str(xi))
            yy.append(float(yi))
            if len(xx) == nrows:
               bool_xy = False
               xydata.append( (xx,yy,cleanstr(pformat),cleanstr(plabel)) )
        if line.startswith("hline ") :
            yvalue, pformat, plabel = line.split()[1:]
            xydata.append( (None,float(yvalue),cleanstr(pformat),cleanstr(plabel)) )
        if line.startswith("vline ") :
            xvalue, pformat, plabel = line.split()[1:]
            xydata.append( (float(xvalue),None,cleanstr(pformat),cleanstr(plabel)) )
        if line.startswith("data ") :
            head, nrows = line.split()[0:2]
            pformat = line.split("'")[1]
            plabel  = line.split("'")[3]
            nrows   = int(nrows)
            xx , yy = [] , []
            if nrows == 0: bool_xy = False
            else         : bool_xy = True
    data = (ptype,xydata,cleanstr(title),cleanstr(xlabel),cleanstr(ylabel))
    return mlabel, data
#--------------------------------------------------#
def read_plotfile(plotfile):
    plotdata = {}
    blocks   = start_end_blocks(plotfile)
    for block in blocks:
        block_lines  = block.split("\n")
        try   : mlabel, data = data_in_block(block_lines)
        except: pass
        plotdata[mlabel] = data
    return plotdata
#==================================================#


#==================================================#
def write_plotfile(plotfile,plotdata):
    # read plot file if exists
    if os.path.exists(plotfile):
        plotdata_file = read_plotfile(plotfile)
        plotdata_file.update(plotdata)
        plotdata = plotdata_file
    # generate string
    string = ""
    for mlabel,(ptype,xydata,title,xlabel,ylabel) in plotdata.items():
        string += "start_%s %s\n"%(ptype,mlabel)
        string += "   title ='%s'\n"%title
        string += "   xlabel='%s'\n"%xlabel
        string += "   ylabel='%s'\n"%ylabel
        for xx,yy,pformat,plabel in xydata:
            # horizontal line
            if   xx is None: string += "   hline %14.5E  '%s'  '%s'\n"%(yy,pformat,plabel)
            # vertical line
            elif yy is None: string += "   vline %14.5E  '%s'  '%s'\n"%(xx,pformat,plabel)
            # points
            else:
               string += "   data  %i  '%s'  '%s'\n"%(len(xx),pformat,plabel)
               for xi,yi in zip(xx,yy):
                   if ptype == "plot": string += "         %14.5E   %14.5E\n"%(xi,yi)
                   if ptype == "bar" : string += "         %s       %14.5E\n"%(xi,yi)
        string += "end_plot\n\n"
    # (re)write plot file
    with open(plotfile,'w') as ff: ff.write(string)
#--------------------------------------------------#
def manage_data_for_plot_weights(ctc,ltemp,dchi):
    plotdata = {}

    itcs = sorted(dchi.keys())
    npts = len(itcs)
    for idx,T in enumerate(ltemp):
        yy     = [dchi[itc][idx] for itc in itcs]
        data   = []
        ptype  = "bar"
        mlabel = "%s_weights_T%0007.2fK"%(ctc,T)
        title  = "Conformational weights for %s at %7.2f K"%(ctc,T)
        xlabel = "conformer"
        ylabel = "weight"
        data.append( (itcs,yy,'','') )
        plotdata[mlabel] = (ptype,data,title,xlabel,ylabel)
    return plotdata
#--------------------------------------------------#
def manage_data_for_plot_mep(tsname,drst,Eref,VadiSpline):

    plotdata = {}

    # data for MEP plot
    data    = []
    mlabel  = "%s_mep"%tsname
    points  = sd.sorted_points(drst,hess=True)
    xx      = [drst[point][0] for point in points]
    yy      = [drst[point][1] for point in points]
    yy      = [(yi-Eref)*pc.KCALMOL for yi in yy]
    title   = "Minimum energy path"
    xlabel  = "reaction coordinate s [bohr]"
    ylabel  = "relative energy [kcal/mol]"
    data.append( (xx,yy,'k--o','low-level') )
    plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)

    # data for Vadi plot
    data    = []
    mlabel  = "%s_vadi"%tsname
    xx      = VadiSpline.xx()
    yy      = [yi*pc.KCALMOL for yi in VadiSpline.yy()]
    xlabel  = "reaction coordinate s [bohr]"
    ylabel  = "relative energy [kcal/mol]"
    title   = "Adiabatic Potential"
    data.append( (xx,yy,'k--o','') )
    plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)

    return plotdata
#--------------------------------------------------#
def manage_data_for_plot_cvt(tsname,ltemp,xx,tgibbs,lscvt,lgamma):

    gibbs,lnew = tgibbs
    plotdata = {}

    # data for Gibbs plot
    for idx,T in enumerate(ltemp):
        data     = []
        mlabel   = "%s_gibbs_T%0007.2fK"%(tsname,T)
        scvt     = lscvt[idx]
        G        = [float(gi*pc.KCALMOL) for gi in gibbs[:,idx]]
        exx, eyy = lnew[idx] # data for -- line
        eyy      = [gi*pc.KCALMOL for gi in eyy]
        xlabel   = "reaction coordinate s [bohr]"
        ylabel   = "relative energy [kcal/mol]"
        title    = "Relative free energy profile at %.2f K"%T
        data.append( (scvt,None,'r--' ,'') ) 
        data.append( (  xx,   G,'ko','') )
        data.append( ( exx, eyy,'--','') )
        plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)

    # data for Gamma_CVT plot
    data    = []
    mlabel  = "%s_cvt"%tsname
    xlabel  = "Temperature [K]"
    ylabel  = r"$\Gamma_{\rm{CVT}}$"
    title   = "CVT variational coefficient"
    data.append( (ltemp,lgamma,'k--o','') )
    plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)

    return plotdata
#--------------------------------------------------#
def manage_data_for_plot_sct(tsname,lkappaZCT,lkappaSCT,svals,mueff,ltemp,\
                  intZCT,intSCT,rteZCT,rteSCT,E0,VAG):
    plotdata = {}

    # data for mu_eff plot
    data    = []
    mlabel  = "%s_mueff"%tsname
    xlabel  = "reaction coordinate s [bohr]"
    ylabel  = "mass [amu]"
    title   = "Effective tunneling mass"
    data.append( (svals,[m*pc.AMU for m in mueff],'k--o','') )
    plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)

    # data for ZCT integrand plot
    for idx,T in enumerate(ltemp):
        data    = []
        mlabel  = "%s_zct_integrand_T%0007.2fK"%(tsname,T)
        xlabel  = "Tunneling Energy [kcal/mol]"
        ylabel  = "ZCT integrand"
        title   = "ZCT integrand at %.2f Kelvin"%T
        xx,yy   = intZCT[idx]
        data.append( (E0*pc.KCALMOL,None,'k--','') )
        data.append( (rteZCT[idx]*pc.KCALMOL,None,'r--','') )
        data.append( (VAG*pc.KCALMOL,None,'k--','') )
        data.append( ([E*pc.KCALMOL for E in xx],yy,'k--o','') )
        plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)

    # data for SCT integrand plot
    for idx,T in enumerate(ltemp):
        data    = []
        mlabel  = "%s_sct_integrand_T%0007.2fK"%(tsname,T)
        xlabel  = "Tunneling Energy [kcal/mol]"
        ylabel  = "SCT integrand"
        title   = "SCT integrand at %.2f Kelvin"%T
        xx,yy   = intSCT[idx]
        data.append( (E0*pc.KCALMOL,None,'k--','') )
        data.append( (rteSCT[idx]*pc.KCALMOL,None,'r--','') )
        data.append( (VAG*pc.KCALMOL,None,'k--','') )
        data.append( ([E*pc.KCALMOL for E in xx],yy,'k--o','') )
        plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)

    # data for ZCT/SCT correction factor plot
    data = []
    mlabel = "%s_zct_sct"%tsname
    xlabel = "Temperature [K]"
    ylabel = r"$\ln(\kappa_{\rm{SCT}})$"
    title  = "SCT transmission coefficient"
    data.append((ltemp,[np.log(yi) for yi in lkappaZCT],'r--D',"ZCT"))
    data.append((ltemp,[np.log(yi) for yi in lkappaSCT],'k--o',"SCT"))
    plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)

    return plotdata
#--------------------------------------------------#
def manage_data_for_plot_rcons(rcname,direction,ltemp,d4plot,units):
    plotdata = {}

    xx1   = [1.0/T for T in ltemp]

    # temperatures for interpolation
    minT   = float(min(ltemp))
    maxT   = float(max(ltemp))
    nT     = 5*len(ltemp)
    dT     = (maxT-minT)/(nT-1)
    ltemp2 = [minT+idx*dT for idx in range(nT)]
    xx2    = [1.0/T for T in ltemp2]
    # data for mu_eff plot
    for rctype,(k,dfit) in d4plot.items():
        data    = []
        mlabel  = "%s_%s-%s"%(rcname,direction,rctype)
        xlabel  = "$1/T$ (K$^{-1}$)"
        ylabel  = "$ln(k)$ ($k$ in $%s$)"%units
        title   = "%s Rate Constant for %s - %s"%(rctype.upper(),rcname,direction)
        # calculated
        yy1 = [np.log(ki) for ki in k]
        data.append( (xx1,yy1,'ko','calculated') )
        #interpolated
        for anatype in dfit.keys():
            coefs, r2 =  dfit[anatype]
            if   anatype == 1: yy2 = [anarc1(T,*coefs) for T in ltemp2]
            elif anatype == 2: yy2 = [anarc2(T,*coefs) for T in ltemp2]
            elif anatype == 3: yy2 = [anarc3(T,*coefs) for T in ltemp2]
            elif anatype == 4: yy2 = [anarc4(T,*coefs) for T in ltemp2]
            else: continue
            # get log
            yy2 = [np.log(ki) for ki in yy2]
            data.append( (xx2,yy2,'--','analytic (%i)'%anatype) )
        # add to dict
        plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)
    return plotdata
#--------------------------------------------------#
def manage_data_for_plot_kmc(data4pdf,fratios,volume,timeunits,logx=True):
    plotdata = {}
    # volume in L
    volume *= pc.ML/1000.0

    # prepare plots
    for stemp, xvalues, yvalues in data4pdf:

        # decide units for time (x-axis)
        tmax = xvalues[-1]*pc.SECOND
        if timeunits == "hr" : tunit = 1./3600*pc.SECOND
        if timeunits == "min": tunit = 1./60  *pc.SECOND
        if timeunits == "s"  : tunit = 1e0    *pc.SECOND
        if timeunits == "ms" : tunit = 1e3    *pc.SECOND
        if timeunits == "mcs": tunit = 1e6    *pc.SECOND
        if timeunits == "fs" : tunit = 1e9    *pc.SECOND
        if timeunits == "ps" : tunit = 1e12   *pc.SECOND
        # add texts
        data = []
        mlabel = "kmc_%sK"%stemp
        if logx: xlabel = "log(time [%s])"%timeunits
        else   : xlabel = "time [%s]"%timeunits
        ylabel = "populations"
        title  = "KMC simulation at T = %s K"%stemp
        # omit first value for logaritmic scale
        if logx: xvalues = [np.log(time*tunit) for time in xvalues[1:]]
        # plot pop of each molecule
        for specie,pops in yvalues.items():
            #concs = [pop/pc.NA/volume for pop in pops]
            if logx: pops = pops[1:]
            data.append( (xvalues,pops,'',specie) )
        plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)

    # plot final ratios
    data    = []
    mlabel  = "kmc_finalratios"
    xlabel  = "T (K)"
    ylabel  = "POP/POP0"
    title   = "KMC simulation - final ratios"
    # plot final pop of each molecule
    ltemp   = sorted(fratios.keys())
    species = sorted(fratios[ltemp[0]].keys())
    xx      = [float(T) for T in ltemp]
    for specie in species:
        yy = [float(fratios[T][specie]) for T in ltemp]
        data.append( (xx,yy,'',specie) )
    plotdata[mlabel] = ("plot",data,title,xlabel,ylabel)
    return plotdata
#--------------------------------------------------#
def plot_mep_dlevel(svals,Ell,Ehl):
    if not PLOT: return
    ll0 = yyll[0]
    hl0 = yyhl[0]
    yyll = [(yi-ll0)*pc.KCALMOL for yi in Ell]
    yyhl = [(yi-hl0)*pc.KCALMOL for yi in Ehl]
    plt.plot(svals,yyll,'-')
    plt.plot(svals,yyhl,'--')
    plt.show()
#==================================================#


