#!/usr/bin/python2.7
'''
*----------------------------------*
| Package    :  common             |
| Module     :  interpolate        |
| Last Update:  2018/10/18 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*


''' 

#==============================================#
from   common.Spline  import Spline
#==============================================#


#=============================================#
#     Functions related to interpolations     #
#=============================================#
def interpolate(xvalues,yvalues,x,d=0):
    spl = Spline(xvalues,yvalues,tension=0.0)
    if d == 0: return spl(x)
    # get derivative
    return spl(x), spl.derivative(x)
#---------------------------------------------#
def interpolate_nones(xx,yy,mode="cubic"):
    '''
    element to correct in yy is given as None
    '''
    npts = len(yy)
    # correct Nones with spline
    if mode == "cubic":
       # create spline
       copy_xx = [x for idx,x in enumerate(xx) if yy[idx] is not None]
       copy_yy = [y for idx,y in enumerate(yy) if yy[idx] is not None]
       spl     = Spline(copy_xx,copy_yy,tension=0.0)
       # interpolate
       for idx in range(npts):
           x,y = xx[idx],yy[idx]
           if y is not None: continue
           yy[idx] = spl(x)
    # correct Nones with linear interpolation
    if mode == "linear":
       for ii in range(npts):
           x,y = xx[ii], yy[ii]
           if y is not None: continue
           idx0 = ii -1
           idx2 = None
           for jj in range(ii,npts):
               if yy[jj] is not None: idx2 = jj; break
           if idx2 is None: continue
           # interpolate
           x0 , y0 = xx[idx0], yy[idx0]
           x2 , y2 = xx[idx2], yy[idx2]
           tangent = (y2-y0)/(x2-x0)
           yy[ii]  = y0 + tangent * (x-x0)
    return yy
#=============================================#

