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

*----------------------------------*
| Package    :  common             |
| Module     :  Logger             |
| Last Update:  2019/04/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the Molecule class
'''

#--------------------------------------------------#
import os
import numpy    as np
import fncs     as fncs
import partfns  as pf
from   dicts    import dpt_im
from   files    import read_gtsfile
from   physcons import AMU, KCALMOL, EV, ANGSTROM, H2CM
from   pgs      import get_pgs
import internal as     intl
import common.Exceptions as Exc
#--------------------------------------------------#

class Molecule():
      '''
      to add: mu
      '''

      def __init__(self):
          self._mform = "empty"

      def __str__(self):
          return self._mform

      def set(self,xcc=None,atonums=None,ch=None,mtp=None,E=None,\
                   gcc=None,Fcc=None,masses=None,pgroup=None,rotsigma=None):
          self._xcc      = xcc
          self._atnums   = atonums
          self._ch       = int(ch)
          self._mtp      = int(mtp)
          self._V0       = E
          self._gcc      = gcc
          self._masses   = masses
          self._pgroup   = pgroup
          self._rotsigma = rotsigma
          # deal with Fcc
          if Fcc is None or Fcc == []:
              self._Fcc = []
          else:
              matrix = np.matrix(Fcc)
              nr, nc = matrix.shape 
              if nr == nc: self._Fcc   = Fcc
              else       : self._Fcc   = fncs.lowt2matrix(Fcc)
          # masses is None?
          if self._masses is None and self._atnums is not None:
             self._masses = fncs.atonums2masses(self._atnums)
          # set ccfreqs to None
          self._ccfreqs  = None
          self._ccFevecs = None

      def prep(self):
          if self._atnums  is not None: self._natoms  = len(self._atnums)
          if self._atnums  is not None: self._symbols = fncs.atonums2symbols(self._atnums)
          if self._symbols is not None: self._mform   = fncs.get_molformula(self._symbols)
          if self._masses  is not None: self._mass    = sum(self._masses)
          if self._mtp     is not None: self._les     = [ (self._mtp,0.0) ] #list of electronic states
          self._fscal = 1.0

      def set_fscal(self,fscal=1.0):
          self._fscal = fscal

      def set_from_gts(self,gtsfile):
          self._gts = gtsfile
          if not os.path.exists(self._gts): return
          xcc,atonums,ch,mtp,E,gcc,Fcc,masses,pgroup,rotsigma,freq_list = read_gtsfile(self._gts)
          # set
          self.set(xcc,atonums,ch,mtp,E,gcc,Fcc,masses,pgroup,rotsigma)
          # some derivated variables
          self.prep()
          # freq list (only for developers)
          if freq_list not in ([],None,()):
             self._ccfreqs = freq_list

      def mod_masses(self,masses):
          self._masses = list(masses)
          self._mass   = sum(self._masses)
          self._pgroup,self._rotsigma = get_pgs(self._atnums,self._masses,self._xcc)

      def apply_imods(self,imods,imasses):
          '''
          example: imods   = ["H2(4,5)","C13(all_C)"]
                   imasses = {"H2":2.0141/AMU, "C13":13.0034/AMU}
          '''
          if imods is None: return
          for imod in imods:
              isymbol = imod.split("(")[0]
              if   isymbol in imasses.keys(): imass = imasses[isymbol]
              elif isymbol in  dpt_im.keys(): imass =  dpt_im[isymbol]
              else:
                 exception = Exc.WrongInIsomass
                 exception._var = isymbol
                 raise exception
              atoms   = imod.split("(")[1].split(")")[0]
              if "all_" in atoms:
                 atype = atoms.split("all_")[1].strip()
                 for idx,symbol in enumerate(self._symbols):
                     if symbol == atype: self._masses[idx] = imass
              else:
                 atoms = [int(ii)-1 for ii in atoms.split(",")]
                 for idx in atoms: self._masses[idx] = imass
          # total mass and point group
          self._mass = sum(self._masses)
          self._pgroup,self._rotsigma = get_pgs(self._atnums,self._masses,self._xcc)

      def mod_les(self,les):
          self._les = les

      def setup(self,mu=1.0/AMU,projgrad=False):
          self._mu = mu
          # shift to center of mass
          self._xcc = fncs.shift2com(self._xcc,self._masses)
          # symmetry
          if self._rotsigma is None:
             self._pgroup,self._rotsigma = get_pgs(self._atnums,self._masses,self._xcc)
          # Molecular case
          if self._natoms != 1:
             # Calculate inertia
             self._itensor = fncs.get_itensor_matrix(self._xcc,self._masses)
             self._imoms, self._rotTs, self._rtype, self._linear = \
                     fncs.get_itensor_evals(self._itensor)
             if self._linear: self._nvdof = 3*self._natoms - 5
             else           : self._nvdof = 3*self._natoms - 6
             # Generate mass-scaled arrays
             self._xms = fncs.cc2ms_x(self._xcc,self._masses,self._mu)
             self._gms = fncs.cc2ms_g(self._gcc,self._masses,self._mu)
             self._Fms = fncs.cc2ms_F(self._Fcc,self._masses,self._mu)
             # calculate frequencies
             if self._ccfreqs is not None: return
             if projgrad: v0 = self._gms
             else       : v0 = None
             data = fncs.calc_ccfreqs(self._Fcc,self._masses,self._xcc,self._mu,v0=v0)
             self._ccfreqs, self._ccFevals, self._ccFevecs = data
             # Scale frequencies
             self._ccfreqs = fncs.scale_freqs(self._ccfreqs,self._fscal)
          # Atomic case
          else:
              self._nvdof = 0
              self._linear = False
              self._xms = list(self._xcc)
              self._gms = list(self._gcc)
              self._Fms = list(self._Fcc)
              self._ccfreqs, self._ccFevals, self._ccFevecs = [], [], []

      def get_imag_main_dir(self):
          ic, fwsign = intl.ics_idir(self._xcc,self._symbols,\
                       self._masses,self._ccfreqs,self._ccFevecs)
          return ic, fwsign

      def icfreqs(self,ics,bool_pg=False):
          if self._natoms != 1:
             ituple = (self._Fcc,self._masses,self._xcc,self._gcc,ics,bool_pg)
             self._icfreqs, self._icFevals, self._icFevecs = intl.calc_icfreqs(*ituple)
          else:
             self._icfreqs  = []
             self._icFevals = []
             self._icFevecs = []
          # scale frequencies
          self._icfreqs = [freq*self._fscal for freq in self._icfreqs]

      def ana_freqs(self,case="cc"):
          if case == "cc":
             # Keep record of imaginary frequencies
             if self._ccFevecs is not None:
                self._ccimag = [ (frq,self._ccFevecs[idx]) for idx,frq in enumerate(self._ccfreqs) if frq < 0.0]
             else:
                self._ccimag = [ (frq,None)                for idx,frq in enumerate(self._ccfreqs) if frq < 0.0]
             # Calculate zpe
             self._cczpes = [fncs.afreq2zpe(frq) for frq in self._ccfreqs]
             self._cczpe  = sum(self._cczpes)
             self._ccV1   = self._V0 + self._cczpe
          if case == "ic":
             # Keep record of imaginary frequencies
             if self._icFevecs is not None:
                self._icimag = [ (frq,self._icFevecs[idx]) for idx,frq in enumerate(self._icfreqs) if frq < 0.0]
             else:
                self._icimag = [ (frq,None)                for idx,frq in enumerate(self._icfreqs) if frq < 0.0]
             # Calculate zpe
             self._iczpes = [fncs.afreq2zpe(frq) for frq in self._icfreqs]
             self._iczpe  = sum(self._iczpes)
             self._icV1   = self._V0 + self._iczpe

      def clean_freqs(self,case="cc"):
          CUTFRQ  = 0.3 # in cm^-1
          # select case
          if case == "cc": freqs = self._ccfreqs
          else           : freqs = self._icfreqs
          # keep track of those to save
          keep = []
          for idx,freq in enumerate(freqs):
              if abs(fncs.afreq2cm(freq)) < CUTFRQ: continue
              keep.append(idx)
          # keep only those > CUTFRQ
          if case == "cc":
             self._ccfreqs  = [self._ccfreqs[idx]  for idx in keep]
             if self._ccFevals is not None:
                self._ccFevals = [self._ccFevals[idx] for idx in keep]
             if self._ccFevecs is not None:
                self._ccFevecs = [self._ccFevecs[idx] for idx in keep]
          if case == "ic":
             self._icfreqs  = [self._icfreqs[idx]  for idx in keep]
             if self._icFevals is not None:
                self._icFevals = [self._icFevals[idx] for idx in keep]
             if self._icFevecs is not None:
                self._icFevecs = [self._icFevecs[idx] for idx in keep]

      def deal_lowfq(self,lowfq={},case="cc"):
          # for Cartesian Coordinates
          if   case == "cc":
             # frequencies were not projected along MEP
             if   self._nvdof - len(self._ccfreqs) == 0:
                for idx,newfreq in lowfq.items():
                    self._ccfreqs[idx] = max(self._ccfreqs[idx],newfreq)
             # frequencies were projected along MEP
             elif self._nvdof - len(self._ccfreqs) == 1:
                for idx,newfreq in lowfq.items():
                    self._ccfreqs[idx-1] = max(self._ccfreqs[idx-1],newfreq)
          # for Internal Coordinates
          elif case == "ic":
             # frequencies were not projected along MEP
             if   self._nvdof - len(self._icfreqs) == 0:
                for idx,newfreq in lowfq.items():
                    self._icfreqs[idx] = max(self._icfreqs[idx],newfreq)
             # frequencies were projected along MEP
             elif self._nvdof - len(self._icfreqs) == 1:
                for idx,newfreq in lowfq.items():
                    self._icfreqs[idx-1] = max(self._icfreqs[idx-1],newfreq)

      def calc_pfns(self,temps,case="cc",fmode=0):
          '''
          fmode = -1 or 0 (0 is default)
          '''
          # Calculate translational partition function (per unit volume)
          ph_tra = np.array([pf.pf_partinbox(self._mass,T) for T in temps])
          # Calculate rotational partition function (Rigid-Rotor)
          if self._natoms > 1:
             pf_rot = np.array([pf.pf_rigidrotor(self._imoms,T,self._rotsigma) for T in temps])
          else:
             pf_rot = np.array([1.0 for T in temps])
          # Calculate vibrational partition function (Harmonic-Oscillator)
          if self._nvdof != 0:
             # remove freq if required
             nf     = self._nvdof + fmode
             if case == "cc": afreqs = list(self._ccfreqs)
             if case == "ic": afreqs = list(self._icfreqs)
             while len(afreqs) > nf: afreqs = afreqs[1:]
             pf_vib = np.array([pf.pf_harmosc(afreqs,T,imag=1E10) for T in temps])
          else:
             pf_vib = np.array([1.0 for T in temps])
          # Calculate electronic partition function
          pf_ele = np.array([pf.pf_electr(self._les,T) for T in temps])
          # Total partition function
          qtot = ph_tra * pf_rot * pf_vib * pf_ele
          if case == "cc": return qtot, self._ccV1, (ph_tra,pf_rot,pf_vib,pf_ele)
          if case == "ic": return qtot, self._icV1, (ph_tra,pf_rot,pf_vib,pf_ele)

      def info_string(self,ib=0):
          root_mass = sum(fncs.symbols2masses(self._symbols))
          string  = "mol. formula      : %s\n"%self._mform
          string += "num atoms         : %i\n"%self._natoms
          string += "num vib dof       : %i\n"%self._nvdof
          string += "charge            : %i\n"%self._ch
          string += "multiplicity      : %i\n"%self._mtp
          string += "total energy (V0) : %.8f hartree\n"%self._V0
          string += "total mass [root] : %.4f amu\n"%(root_mass *AMU)
          string += "total mass        : %.4f amu\n"%(self._mass*AMU)
          if self._pgroup   is not None: string += "point group       : %s\n"%(self._pgroup)
          if self._rotsigma is not None: string += "rot sym num       : %i\n"%(self._rotsigma)
          string += "Cartesian coordinates (Angstrom):\n"
          for at,symbol in enumerate(self._symbols):
              mass   = self._masses[at]*AMU
              x,y,z  = fncs.xyz(self._xcc,at)
              x *= ANGSTROM
              y *= ANGSTROM
              z *= ANGSTROM
              string += "  %2s   %+10.6f  %+10.6f  %+10.6f  [%7.3f amu]\n"%(symbol,x,y,z,mass)

          #if True:
          try:
              str2  = "moments and product of inertia (au):\n"
              if len(self._imoms) == 1:
                 str2 += "        %+10.3E\n"%self._imoms[0]
              if len(self._imoms) == 3:
                 prodinert = self._imoms[0]*self._imoms[1]*self._imoms[2]
                 dataline = (self._imoms[0],self._imoms[1],self._imoms[2],prodinert)
                 str2 += "        %+10.3E  %+10.3E  %+10.3E  [%10.3E]\n"%dataline
              string += str2
          except: pass

          try:
              str2  = "vibrational frequencies [1/cm] (scaled by %.3f):\n"%self._fscal
              for idx in range(0,len(self._ccfreqs),6):
                  str2 += "  %s\n"%("  ".join("%8.2f"%fncs.afreq2cm(freq) \
                                      for freq in self._ccfreqs[idx:idx+6]))
              if len(self._ccfreqs) != 0: string += str2
          except: pass

          try:
              str2  = "individual zero-point energies [kcal/mol]:\n"
              for idx in range(0,len(self._cczpes),6):
                  str2 += "  %s\n"%("  ".join("%8.2f"%(zpe*KCALMOL) \
                                      for zpe in self._cczpes[idx:idx+6]))
              zpe_au   = self._cczpe
              zpe_kcal = self._cczpe * KCALMOL
              zpe_eV   = self._cczpe * EV
              zpe_cm   = self._cczpe * H2CM
              str2 += "zero-point energy (zpe) : %+14.8f hartree  = \n"%zpe_au
              str2 += "                          %+14.2f kcal/mol = \n"%zpe_kcal
              str2 += "                          %+14.2f eV       = \n"%zpe_eV
              str2 += "                          %+14.2f cm^-1 \n"%zpe_cm
              #str2 += "zero-point energy (zpe) : %+14.8f hartree = %.2f kcal/mol = %.3f eV\n"%(zpe_au,zpe_kcal,zpe_eV)
              str2 += "total energy + zpe (V1) : %+14.8f hartree\n"%self._ccV1
              if self._cczpe != 0.0: string += str2
          except: pass

          # add blank spaces
          string = "\n".join([" "*ib+line for line in string.split("\n")])
          return string

