#!/usr/bin/python2.7
'''
*----------------------------------*
| Package    :  common             |
| Module     :  read_gauout        |
| Last Update:  2019/01/26 (Y/M/D) |
| Main Author:  Antonio Fernandez  |
*----------------------------------*

'''

#=======================================================#
import numpy           as     np
#-------------------------------------------------------#
import Exceptions      as     Exc
#-------------------------------------------------------#
from   common.pgs      import get_pgs
from   common.fncs     import flatten_llist
from   common.fncs     import atonums2masses
from   common.physcons import ANGSTROM
#=======================================================#

def archive(filename,natoms,start_archive,end_archive,\
            pos_beginfile,pos_endfile,pos_nosym,pos_geom1,pos_geom2,pos_zmat1): 
   # lists that returns:
   #  - list_level: a list of all levels of calculation
   #  - list_energy: the energy of each of the levels
   #  - ccgeom: the geometry in cartesian coordinates
   #            1) From the output file: cartesian 
   #                OR if Lxyz
   #            2) from the archive: xyz
   #  - atomic_number: a list with the atomic numbers
   #  - list_zmat: (if Lzmat) a list with the z-matrix updated 
   #  - hess_full: (if Lhess) a list with the full hessian 
   #  - grad     :   ""       a list with the gradient   

   ndim=3*natoms

   bar='\\'
   bar2='\\\\'
   fin='@'
   ver='Version'
   str_en1='State='
   str_en2='RMSD='
   hag='#'
   str_imag='NImag='
   str_hessian='Frequencies'
   Lzmat=False
   Lhess=False
   Lxyz=False
   lista_archive=[]
   # The program always reads the Z-matrix or Input orientation
   #if pos_nosym > pos_beginfile: pos_geom=pos_geom1
   #else: pos_geom=pos_geom2
   pos_geom = pos_geom1

   # Z-matrix?
   if pos_zmat1 > 0: 
      Lzmat=True
      pos_zmat2=start_archive
      
   # Reads the last archive file and removes returns
   # Also reads the geometry and Z-matrix if available 
   # It removes the dummy atoms from the geometry
   nlin=0
   idummy = 0
   list_geom=[]
   list_zmat=[]
   atomic_number=[]
   atomic_symbol=[]
   cartesian=[]
   tmp_list=[]
   filex=open(filename,'r') 
   for linea in filex:
      list_dummy=[]
      nlin+=1
      if nlin >= pos_beginfile and nlin <= pos_endfile:
         if str_hessian in linea: Lhess = True
      linea=linea.strip()
      if nlin >= pos_geom+5 and nlin <= pos_geom+4+natoms+idummy:
         list_dummy = linea.split()
         if int(list_dummy[1]) > 0:
            list_geom.extend(linea.split())
         else:
            idummy += 1
      if nlin >= start_archive and nlin <= end_archive:
         lista_archive.append(linea)
      if Lzmat and nlin > pos_zmat1 and nlin < pos_zmat2:
	    list_zmat.append(linea)
   filex.close()
   
   ncol_geom=6
   for i in range(natoms):
      atomic_number.append(int(list_geom[i*ncol_geom+1]))
      tmp_list=[float(list_geom[i*ncol_geom+3]),float(list_geom[i*ncol_geom+4]),float(list_geom[i*ncol_geom+5])]
      cartesian.append(tmp_list) 


   # Merges the list into one element
   length_lista=len(lista_archive)
   lista_archive[0:length_lista]=[''.join(lista_archive[0:length_lista])]
   length_lista=len(lista_archive)
   str_lista=lista_archive[0]

   # Command line   
   pcom1=str_lista.find(hag,0)
   pcom2=str_lista.find(bar2,pcom1)
   commands=str_lista[pcom1:pcom2]

   # Charge and multiplicity
   pchm1=str_lista.find(bar2,pcom2+len(bar2))
   pchm2=str_lista.find(bar,pchm1+len(bar2))
   t2=str_lista[pchm1+len(bar2):pchm2]
   chmul=map(int, t2.split(',') )
   charge=chmul[0]
   mult=chmul[1]

   # zmatrix or cartesian coordinates
   pgm1=pchm2+1
   pgm2=str_lista.find(ver,0)
   str_geom=str_lista[pgm1:pgm2]
   geom=map( str, str_geom.split('\\') )
   if len(geom[0]) <= 4:
      Lzmat=True
      list_zmat=geom
   else:
      Lxyz=True
      xyz_tmp=[]
      tmp_list=[]
      xyz=[]
      ncol_xyz=4
      for i in range(natoms):
	 xyz_tmp.extend(geom[i].split(','))
      for i in range(natoms):
	 atomic_symbol.append(xyz_tmp[ncol_xyz*i])
	 tmp_list=[float(xyz_tmp[i*ncol_xyz+1]),float(xyz_tmp[i*ncol_xyz+2]),float(xyz_tmp[i*ncol_xyz+3])]
	 xyz.append(tmp_list) 
   if Lxyz:
      ccgeom=xyz
   else:
      ccgeom=cartesian 


   # Energies & other info
   list_level=[]
   list_energy=[]
   #kk=str_lista.find(str_en1)
   #if kk < 0: str_first = pgm2
   #else: str_first = str_en1
   pener1x=str_lista.find(ver,0)
   pener1y=str_lista.find(str_en1,0)
   pener1 = max(pener1x,pener1y)
   pener2=str_lista.find(bar,pener1)
   pener3=str_lista.find(str_en2,pener1)
   t3=str_lista[pener2+1:pener3-1]
   t4=t3.replace(bar,' ').replace('=',' ').split()
   num_energies=0
   for i in range(0,len(t4),2):
      # avoids the incorporation of S**2 in the list of energies
      # for openshell systems
      t4sub = t4[i]
      if t4sub[0] != 'S':   
         num_energies +=1
         list_level.append(t4[i])
         list_energy.append(float(t4[i+1]))

   #Hessian
   if Lhess:
      phes1=str_lista.find(str_imag,0)
      phes2=str_lista.find(bar2,phes1)
      phes3=str_lista.find(bar2,phes2+len(bar2))
      num_imag=int(str_lista[phes1+len(str_imag):phes2])
      str_hess=str_lista[phes2+len(bar2):phes3]
      hess_data=map( float, str_hess.split(',') )
      hess_lt = np.zeros((ndim,ndim))
      indices = np.tril_indices(ndim)
      hess_lt[indices] = hess_data
      hess_full=hess_lt
      # Full hessian
      for i in range(ndim):
       for j in range(i+1):
	   hess_full[j][i] = hess_full[i][j]
      
   else:
      num_imag=-1
      pass

   #Gradient
   if Lhess > 0:
      pgrd1=phes3+len(bar2)
      pgrd2=str_lista.find(fin,0)
      str_grd=str_lista[pgrd1:pgrd2-3]
      grad=map( float, str_grd.split(',') )
   
   if not Lzmat: list_zmat=[0] 
   if not Lhess:
      hess_data=[0]
      hess_full=[0]
      grad=[0]

   return num_imag,charge,mult,commands,list_level,list_energy,atomic_number,ccgeom,\
          list_zmat,hess_data,hess_full,grad


def read_gauout(filename):

    # check extension
    end_out = (filename.lower()).endswith(".out")
    end_log = (filename.lower()).endswith(".log")
    if (not end_out) and (not end_log):
        raise Exc.FileType(Exception)

    key='GINC'
    fin='@'
    str_end='Normal termination'
    str_atoms='NAtoms'
    str_geom1a='Z-Matrix orientation'
    str_geom1b='Input orientation'
    str_geom2='Standard orientation'
    str_nosym='Nosym'
    str_zmat='Final structure in terms of initial Z-matrix'
    start_archives=[]
    end_archives=[]
    str_endfiles=[]
    latom=[]
    idx_last_zmat = -1
    list_zmat_backup = []

    # More then one archive file for instance Link1
    filex=open(filename,'r') 
    nlin=0
    for linea in filex:
       nlin+=1
       if str_end in linea: str_endfiles.append(nlin)
       if key in linea: start_archives.append(nlin)
       if fin in linea: end_archives.append(nlin)
       if str_atoms in linea : 
          latom=linea.split() 
    filex.close()
    # It does not consider an archive file right after the main archive file
    # For some reason gaussian sometimes print one archive after the other
    str_tmp = str_endfiles[:] 
    num_tmp = len(str_tmp)
    if num_tmp > 1:
       for i in range(num_tmp-1,0,-1):
          if abs(str_tmp[i-1]-str_tmp[i]) < 300:
             del str_endfiles[i]
             del start_archives[i]
             del end_archives[i]
    #
    #
    number_archives=len(start_archives)
    if number_archives == 0: return -1
    str_beginfiles=str_endfiles[:]
    str_beginfiles.insert(0,0)
    str_beginfiles.pop()

    natoms = int(latom[1])
    if number_archives == 1 and start_archives[0] == 0: return

    str_geom1s=[]
    str_geom2s=[]
    str_nosyms=[]
    str_zmats=[]
    for i in range(number_archives):
       str_geom1s.append(0)
       str_geom2s.append(0)
       str_nosyms.append(0)
       str_zmats.append(0)
       start_archive=start_archives[i]
       end_archive=end_archives[i]
       pos_beginfile=str_beginfiles[i]
       pos_endfile=str_endfiles[i]
       filex=open(filename,'r') 
       nlin=0
       for linea in filex:
          nlin+=1
          if nlin >= pos_beginfile and nlin <= pos_endfile :
             if str_nosym.lower() in linea.lower(): str_nosyms[i]=nlin
             if str_geom1a in linea or str_geom1b in linea : str_geom1s[i]=nlin
             if str_geom2 in linea : str_geom2s[i]=nlin
             if str_zmat in linea : str_zmats[i]=nlin
       filex.close()

       # position of the geometries
       pos_nosym=str_nosyms[i]
       pos_geom1=str_geom1s[i]
       pos_geom2=str_geom2s[i]
       pos_zmat1=str_zmats[i]
       num_imag,charge,mult,commands,list_level,list_energy,atomic_number,ccgeom,list_zmat,\
           hess_data,hess_full,grad=\
           archive(filename,natoms,start_archive,end_archive,\
           pos_beginfile,pos_endfile,pos_nosym,pos_geom1,pos_geom2,pos_zmat1)

       # Saves the last available z-matrix
       if str_zmats[i] != 0:
          idx_last_zmat = i
          list_zmat_backup = list_zmat[:]

    # atomic mass
    atomic_mass = atonums2masses(atomic_number)

    ccgeom = flatten_llist(ccgeom)

    energy = list_energy[-1]
    level  = list_level[-1]
    # from Angstrom to Bohr
    ccgeom = [xx/ANGSTROM for xx in ccgeom]

    return ccgeom, atomic_number, charge, mult, energy, grad, hess_data, atomic_mass, level


