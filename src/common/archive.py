
import numpy as np
import sys
import os

def get_pgs(atom_num,atom_mass,geom_xyz,toldist=0.05,tolsym=3e-2,epsilon=3e-2):
  '''
  This module finds the symmetry point
  group symmetry of a molecule 
  
  *---------------------------------------*
  | Main Author:  Antonio Fernandez-Ramos |
  | Last Update:  May 10th 2017 (by DFC)  |
  *---------------------------------------*
  '''
  
  #  -------Parameters----------
  # atomic numbers
  # atomic masses
  # Geometry as [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3] ....
  # toldist : tolerance of difference of distances
  # tolsym : tolerance in the symmetry operation
  # epsilon : tolerance in the cross products
  #           usually used to check if two C2 axis are parallel

  # /////  Matrix representation of symmetry operators \\\\\
  
  # in case geom_xyz is a 3N list 
  if len(geom_xyz) == 3*len(atom_num):
     xyz = []
     for idx in range(0,len(geom_xyz),3):
         x,y,z = geom_xyz[idx+0:idx+3]
         xyz.append( [x,y,z] )
     geom_xyz = xyz

  natom = len(atom_num)

  if natom == 1: return 'C1',1

  def inversion():
      i_center = ([[-1.,0.,0.],
                   [0.,-1.,0.],
                   [0.,0.,-1.]])
      return i_center
  
  def reflex(uvw,main_plane):
  # reflection about a plane perpendicular to uvw
    plane = np.zeros( (3,3) )
    if main_plane == 'xy': uvw = [0.,0.,1.]
    if main_plane == 'xz': uvw = [0.,1.,0.]
    if main_plane == 'yz': uvw = [1.,0.,0.]
    x = uvw[0]
    y = uvw[1]
    z = uvw[2]
    plane[0,0] = 1.-2.*x**2 
    plane[1,1] = 1.-2.*y**2 
    plane[2,2] = 1.-2.*z**2 
    plane[0,1] = -2.*x*y
    plane[0,2] = -2.*x*z
    plane[1,2] = -2.*y*z
    plane[1,0] = plane[0,1] 
    plane[2,0] = plane[0,2] 
    plane[2,1] = plane[1,2] 
    return plane
  
  def Cngen(n,u):
  # rotacion about the u vector
     Cngen = np.zeros( (3,3) )
     t=2.*np.pi/float(n)
     ux = u[0]
     uy = u[1]
     uz = u[2]
     cosu = np.cos(t)
     sinu = np.sin(t)
     Cngen[0,0] = cosu + ux**2*(1.-cosu)
     Cngen[1,1] = cosu + uy**2*(1.-cosu)
     Cngen[2,2] = cosu + uz**2*(1.-cosu)
     Cngen[0,1] = ux*uy*(1.-cosu)-uz*sinu
     Cngen[0,2] = ux*uz*(1.-cosu)+uy*sinu
     Cngen[1,2] = uy*uz*(1.-cosu)-ux*sinu
     Cngen[1,0] = ux*uy*(1.-cosu)+uz*sinu
     Cngen[2,0] = ux*uz*(1.-cosu)-uy*sinu
     Cngen[2,1] = uy*uz*(1.-cosu)+ux*sinu
     return Cngen

### ///// END \\\\
  
  def cdmass(mass,xyz,natom):
   # center of masses
     cdm = [0.,0.,0.]
     tot_mass = sum(mass)
     for i in range(0,natom):
       cdm[0] += mass[i]*xyz[i,0]/tot_mass
       cdm[1] += mass[i]*xyz[i,1]/tot_mass
       cdm[2] += mass[i]*xyz[i,2]/tot_mass
     for i in range(0,natom):
       xyz[i,0] += - cdm[0]
       xyz[i,1] += - cdm[1]
       xyz[i,2] += - cdm[2]
     return xyz
  
  def I_diff(fx,fy,tol):
      abs_rel = np.fabs(fx-fy)/fy*100.
      if abs_rel <= tol: 
        return True
      else:
        return False
  
  def calc_tensor(mass,xyz,natom):
    # tensor of inertia 
      tensor_iner = np.zeros( (3,3) )
      for i in range(0,natom):
        tensor_iner[0,0] += mass[i]*(xyz[i,1]**2+xyz[i,2]**2)
        tensor_iner[1,1] += mass[i]*(xyz[i,0]**2+xyz[i,2]**2)
        tensor_iner[2,2] += mass[i]*(xyz[i,0]**2+xyz[i,1]**2)
        tensor_iner[0,1] += -mass[i]*xyz[i,0]*xyz[i,1]
        tensor_iner[0,2] += -mass[i]*xyz[i,0]*xyz[i,2]
        tensor_iner[1,2] += -mass[i]*xyz[i,1]*xyz[i,2]
      tensor_iner[1,0] = tensor_iner[0,1]
      tensor_iner[2,0] = tensor_iner[0,2]
      tensor_iner[2,1] = tensor_iner[1,2]
      #eigvals, eigvecs = np.linalg.eigvalsh(tensor_iner)
      Im,Ivec =  np.linalg.eigh(tensor_iner)
      Ivec = Ivec.transpose()
      return Im,Ivec
  
  
  def dist_mat(xyz,natom):
  # Matrix of distances
     dmx = []
     for i in range(0,natom):
       vtmp = xyz[i]
       for j in range (0,natom):
         vtmp2 = xyz[j]
         vx=vtmp-vtmp2
  # Norm + storage
         dmx.append(np.linalg.norm(vx))
     dmx=np.reshape(dmx,(natom,natom))
     return dmx
  
  def dist_c(v1,v2):
     L_novec = False
     vxn = []
     vx = (v1+v2)/2
     vxn=np.linalg.norm(vx)
     if vxn == 0.: 
       L_novec = True
       return vx,L_novec
     vx /= vxn
     return vx,L_novec 
     
  def dist_v(v1,v2):
     L_nosig = False
     vxn = []
     vx = (v1-v2)
     vxn=np.linalg.norm(vx)
     if vxn == 0.: 
       L_nosig = True
       return vx,L_nosig
     vx /= vxn
     return vx,L_nosig
    
  
  def get_sea(dict_sea,mat_geo,idx):
  # Mass and geometry of each SEA
    coor_sea = []
    mass_sea = []
    idum=len(dict_sea[idx])
    for i in range(idum):
          k = dict_sea[idx][i]
          coor_sea.append(mat_geo[k])
          mass_sea.append(atom_mass[k])
    coor_sea=np.reshape(coor_sea,(idum,3))
    return idum,mass_sea,coor_sea
  
  def compare_geom(A,B):
    nrows, ncols = A.shape
    tot_error = 0.0
    indices_comparison = []

    for row_i in range(nrows):
        Arow = A[row_i,:]
        diff    = float("inf")
        idx     = None
        for row_j in range(nrows):
            Brow = B[row_j,:]
            error = np.linalg.norm(Arow - Brow)
            if error < diff:
               diff = error
               idx  = row_j
        tot_error += diff
        indices_comparison.append( (row_i,idx)  )
    tot_error = 1.0 * tot_error / nrows
    return tot_error 

  def get_sym(dict_sea,mat_geo,mat_sym,num_set,strsym,tolsym):
    # returns True or False, i.e. if the molecule is invariant after a given
    # symmetry operation (strsym) 
    total_error = 0.
    mat_tmp = mat_geo.T 
    mat_tmp = mat_sym*mat_tmp
    mat_tmp = mat_tmp.T
    L_sym = False
    for i in range(num_set):
      idx = i 
      isea,mass_sea,coor_sea =  get_sea(dict_sea,mat_geo,idx)
      isea,mass_tmp,coor_tmp =  get_sea(dict_sea,mat_tmp,idx)
      coor_sea = np.array(coor_sea)
      error_sym_oper = compare_geom(coor_sea,coor_tmp)
      total_error += error_sym_oper
    total_error /= num_set
#    print 'Error ',total_error
    if total_error < tolsym:
      L_sym = True

    return L_sym  
  
  def get_c2(dict_sea,mat_geo,num_set,uvw,L_round,epsilon):
    #
    # uvw is for all type of molecules except spheric
    # L_round = True for spherical molecules 
    # search for C2 axis: 1) passing through the middle of a bond
    #                     2) passing through atoms
    coor_sea = []
    mass_sea = []
    new_c2 = 0
    new_c2_red = 0  
    list_tot = []
    list_red = []
    dxc2 = []

    for l in range(num_set):
      idx = l
      isea,mass_sea,coor_sea = get_sea(dict_sea,mat_geo,idx)
      if isea > 1:
        for i in range(isea):
          for j in range(i+1,isea):
            dx,L_novec = dist_c(coor_sea[i],coor_sea[j])
            if L_novec:
              continue
            elif L_round:
              dxm.append(dx)
              new_c2_red += 1
            else:
              vcros = np.cross(dx,uvw) 
              vcrosn=np.linalg.norm(vcros)
              if vcrosn > 1.-epsilon: 
                new_c2_red += 1
                dxm.append(dx)
          vnorm=np.linalg.norm(coor_sea[i])
          if vnorm > epsilon:
            dxm.append(coor_sea[i]/vnorm)
            new_c2_red += 1
  #
  # ---remove the redundancies from the C2 axes
  #
    for i in range(new_c2_red): list_tot.append(i)  
    set_tot = set(list_tot)
    if new_c2_red > 1:
      for i in range(new_c2_red):
        vdxi = dxm[i]
        for j in range(i+1,new_c2_red): 
          vdxj = dxm[j]
          vcros = np.cross(vdxi,vdxj) 
          vcrosn=np.linalg.norm(vcros)
          if j in list_red:
            continue
          else:
            if vcrosn < epsilon: 
              list_red.append(j)
    set_red = set(list_red)
    set_not_red = set_tot - set_red
    dxc2.append([dxm[i] for i in set_not_red])
    new_c2 = len(set_not_red)
    dxc2=np.reshape(dxc2,(new_c2,3))
    return new_c2,dxc2

  def get_sigma_d(dict_sea,mat_geo,num_set,mat_c2_axis,numc2,tolsym):
  #
  #  Diagonal planes: considers that are perpendicular to the
  #    vectors defined by the C2 axes and the coordinate of
  #    an atom of the SEA 
  #    Important for Td molecules
  #
    coor_sea = []
    mass_sea = []
    mat_sigma_d = [] 
    num_sigma_d = 0
    vdxi = []
    vdxj = []
    for l in range(num_set):
      idx = l
      isea,mass_sea,coor_sea = get_sea(dict_sea,mat_geo,idx)
      if isea > 1:
        for i in range(isea):
          vdxi = coor_sea[i]
          for j in range(numc2):
            vdxj = mat_c2_axis[j]
            vnormal = np.cross(vdxi,vdxj)
            vnorm=np.linalg.norm(vnormal)
            if vnorm > tolsym:
              vnormal /= vnorm 
              mat_sigma_d = reflex(vnormal,' ')
              L_sigma_d = get_sym(dict_sea,mat_geo,mat_sigma_d,num_set,'sigma_d',tolsym) 
              if L_sigma_d:
                num_sigma_d += 1
            else:
              continue
    return num_sigma_d 
    
  
  #####initial variables##################
  
  L_cn = False
  L_c2_p = False
  L_c2_sp = False
  L_sigma_v = False
  L_sigma_h = False
  L_sigma_x = False
  L_i = False
  L_s2n = False
  L_cs = False
  cnrot = []
  cn_from_I = [1]
  sigma_h = []
  sigma_v = []
  sigma_x = []
  sn_mat = []
  i_center = inversion()
  udum = []
  uvw = [0.,0.,0.]
  axis = ['x','y','z']
  plane = ['yz','xz','xy']
  dxm = [] 
  
  
  mat_geo = np.matrix(geom_xyz)
  # Puts the molecule in the center of masses
  mat_geo = cdmass(atom_mass,mat_geo,natom)
  # Evaluates the tensor of inertia
  Ipmol,Ipvec = calc_tensor(atom_mass,mat_geo,natom)
  # Rotates the molecule so Ia, Ib and Ic coincide with i,j,k
  mat_geo = Ipvec * mat_geo.T
  mat_geo = mat_geo.T
  
  
  dist = dist_mat(mat_geo,natom)
  distv = []
  # sum the distances of each row
  for i in range(natom):
     rkk = sum(dist[i,:])/float(natom)
     distv.append(rkk)
  
  # now it checks for possible SEA
  # index of SEAs stored in the dictionary
  dict_sea = { }
  idx = -1
  sea_at_idx = []
  for i in range(natom):
     next = False
     if i not in sea_at_idx: 
       idx += 1
     for j in range(i,natom):
       
       if atom_mass[i] == atom_mass[j] and atom_num[i] == atom_num[j] and np.fabs(distv[i] - distv[j])< toldist and j not in sea_at_idx:
         sea_at_idx.append(j)
         dict_sea[idx]  = dict_sea.get(idx,[]) + [j]
  
  # number of SEAs
  # maximum number of atoms in SEA
  num_set = idx+1
  atoms_sea = []
  for i in range(num_set):
    atoms_sea.append(len(dict_sea[i])) 
  atom_sea_max = max(atoms_sea)
  
  
  # Is it linear?
  if np.fabs(Ipmol[0]) <= epsilon:
  # Linear molecule
     linear = True
     L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym)
     if L_i == True:
       return 'Dinfv',2
     else:
       return 'Cinfv',1
  
  # Non linear molecule
      # spheric symmetry (Ia = Ib = Ic)
      # achatado : oblate (Ia = Ib < Ic)
      # alargado : prolate (Ia < Ib = Ic)
      # asymmetric (Ia != Ib != Ic)
  L_sphera = False
  L_asym = False
  L_oblate = I_diff(Ipmol[0],Ipmol[1],epsilon*10)
  L_prolate = I_diff(Ipmol[1],Ipmol[2],epsilon*10)
  if L_oblate and L_prolate:
    L_sphera = True
    L_asym = False
  if L_sphera:
    L_prolate = False
    L_oblate = False
  if not L_oblate and not L_prolate and not L_sphera: 
    L_asym = True

  #
  # --- SPHERICAL MOLECULES
  #
  if L_sphera: 
    # ---C2 axes
    uvw_asym = [0.,0.,0]
    new_c2,dxc2 = get_c2(dict_sea,mat_geo,num_set,uvw_asym,L_sphera,epsilon)
    c2sp = []
    dxc2_save = []
    number_c2 =0
    for i in range(new_c2):
      c2sp = Cngen(2,dxc2[i])
      L_c2_sp=get_sym(dict_sea,mat_geo,c2sp,num_set,'C2',tolsym) 
      if L_c2_sp:
        number_c2 += 1
        dxc2_save.append(dxc2[i]) 
    # ---diagonal planes
    if number_c2 <= 3:
      num_sigma_d = get_sigma_d(dict_sea,mat_geo,num_set,dxc2_save,number_c2,tolsym)
    # ---inversion
    L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym)
    if number_c2 == 3:
      if L_i:
        return 'Th',12
      elif num_sigma_d >= 6: 
        return 'Td',12 
      else:
        return 'T',12 
    if number_c2 == 9:
      if L_i:
        return 'Oh',24
      else:
        return 'O',24 
    if number_c2 > 9:
      if L_i:
        return 'Ih',60
      else:
        return 'I',60 
  
  #
  # --- SYMMETRIC MOLECULES (PROLATE or OBLATE)
  #
  #      Cn axis search
  #
  if L_oblate or L_prolate:
      L_cn = False
      L_cn_tmp = False
      L_cn2_tmp = False
      if L_oblate: 
         main_axis = axis[2]
         main_plane = plane[2]
         uvw = [0.,0.,1.]
      if L_prolate: 
         main_axis = axis[0]
         main_plane = plane[0]
         uvw = [1.,0.,0.]
      for nrot in range(6,1,-1):
        cn_axis = nrot
        cad_axis = 'c'+str(cn_axis)+main_axis
        cnrot =  Cngen(cn_axis,uvw)
        L_cn=get_sym(dict_sea,mat_geo,cnrot,num_set,cad_axis,tolsym) 
        if L_cn: 
          max_cn_axis = cn_axis
          break
      if not L_cn:
        L_asym = True
  
  #
  # --- ASYMMETRIC MOLECULES
  #
  #      C2 axis search
  #
  if L_asym:
    max_cn_axis = 1
    for iaxis in range(3):
      init_axis = axis[iaxis]
      init_plane = plane[iaxis]
      if iaxis == 0: uvw = [1.,0.,0.]
      if iaxis == 1: uvw = [0.,1.,0.]
      if iaxis == 2: uvw = [0.,0.,1.]
      cn_axis = 2
      cad_axis = 'c'+str(cn_axis)+init_axis
      cnrot = []
      cnrot =  Cngen(cn_axis,uvw)
      L_cn_tmp=get_sym(dict_sea,mat_geo,cnrot,num_set,cad_axis,tolsym) 
      if L_cn_tmp:
        L_cn = L_cn_tmp
        max_cn_axis = 2
        main_axis = init_axis
        main_plane = init_plane
        cad_plane = 'sigma_'+main_plane
        dxm.append(uvw) 
        uvw_asym = uvw
  
  # --- C2 axis perpendicular to the main axis with redundancies
  # ... Check for vertical planes too
  # --- If not C2 then looks for a Cs plane

  #
  # ---finds C2 axis 
  #
  if L_cn:
    if L_prolate or L_oblate: uvw_asym = uvw 
    new_c2,dxc2 = get_c2(dict_sea,mat_geo,num_set,uvw_asym,L_sphera,epsilon)
    new_c2_p = 0
    if new_c2 > 0:
      for i in range(new_c2):
        cnrot = Cngen(2,dxc2[i])
        L_kk=get_sym(dict_sea,mat_geo,cnrot,num_set,'C2p',tolsym) 
        if L_kk: new_c2_p += 1
    if new_c2_p >= max_cn_axis : L_c2_p = True
  #
  # ---finds vertical planes sigma_v
  #
    coor_sea = []
    mass_sea = []
    num_sigma_v = 0
    for l in range(num_set):
     idx = l
     isea,mass_sea,coor_sea =  get_sea(dict_sea,mat_geo,idx)
     if isea > 1:
      for i in range(isea):
        for j in range(i+1,isea):
          dv,L_nosig = dist_v(coor_sea[i],coor_sea[j])
          if L_nosig:
            continue
          else: 
            sigma_v = reflex(dv,' ')
            L_sigma_tmp=get_sym(dict_sea,mat_geo,sigma_v,num_set,' ',tolsym) 
            if(L_sigma_tmp): num_sigma_v += 1
    
  
    if num_sigma_v > 0: L_sigma_v = True
  #
  #--- look for horizontal reflection
  #
    sigma_h = np.matrix(reflex(udum,main_plane))
    cad_plane = 'sigma_h'
    L_sigma_h=get_sym(dict_sea,mat_geo,sigma_h,num_set,cad_plane,tolsym) 
  #--- look for inversion
    L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym)
  # ---look for S2n symmetry operation
    sn_axis = 2*max_cn_axis
    cnrot =  Cngen(sn_axis,uvw)
    snrot = sigma_h * cnrot
    cad_saxis = 's'+str(sn_axis)
    L_s2n=get_sym(dict_sea,mat_geo,snrot,num_set,cad_saxis,tolsym) 
  else:
   #-- looking for a reflection plane
    for iaxis in range(3):
      init_axis = axis[iaxis]
      init_plane = plane[iaxis]
      cad_cs = 'Cs in axis ',axis[iaxis]
      csplane = reflex(udum,init_plane)
      L_cs=get_sym(dict_sea,mat_geo,csplane,num_set,cad_cs,tolsym)
      if L_cs:
        break
    uvw = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
    num_sigma_cs = 0
    num_sigma_cs = get_sigma_d(dict_sea,mat_geo,num_set,uvw,3,tolsym)
    if num_sigma_cs > 0: L_cs = True
   #-- looking for inversion
    L_i=get_sym(dict_sea,mat_geo,i_center,num_set,'inversion',tolsym)

# Any Cn axis (n >=2)?
#    YES
#       Any C2 axis perpendicular to the Cn? 
#         YES  
#            Horizontal plane?
#               YES -> Dnh
#               NO
#                  Improper rotation?
#                     YES -> Dnd
#                     NO  -> Dn
#         NO 
#            Horizontal plane?
#               YES -> Cnh
#               NO  
#                  Vertical planes?
#                     YES -> Cnv
#                     NO
#                        Improper rotation?
#                           YES -> S2n
#                           NO  -> Cn
#    NO  
#      Any plane?
#         YES -> Cs
#            Inversion?
#               YES -> Ci
#               NO  -> C1

  if L_cn:
    if L_c2_p:
      if L_sigma_h:
        return 'D'+str(max_cn_axis)+'h',2*max_cn_axis
      else:
        if L_s2n:
          return 'D'+str(max_cn_axis)+'d',2*max_cn_axis
        else:
          return 'D'+str(max_cn_axis),2*max_cn_axis
    else:
      if L_sigma_h:
        return 'C'+str(max_cn_axis)+'h',max_cn_axis
      else: 
        if L_sigma_v:
          return 'C'+str(max_cn_axis)+'v',max_cn_axis
        else:
          if L_s2n:
            return 'S'+str(sn_axis),max_cn_axis
          else:
            return 'C'+str(max_cn_axis),max_cn_axis
  else:
    if L_cs: 
      return 'Cs',1
    else:
      if L_i: 
        return 'Ci',1
      else:
        return 'C1',1

def write_gtsfile(xyz_list , atonum_list , charge, multiplicity, energy, \
                   pgroup, rotsigma , g_list , F_list, gtsfile, freqs=None):

    natoms   = len(atonum_list)
    xyz_list = xvecformat(xyz_list,natoms,out="Nx3")
    if g_list is not None: g_list = xvecformat(g_list,natoms,out="Nx3")


    gts = open(gtsfile,"w")

    # Write atomic numbers and cartesian coordinates
    gts.write("# Atomic number and non-scaled cartesian coordinates [bohr]\n")
    gts.write("start_cc\n")
    for atonum, xyz in zip(atonum_list,xyz_list):
        x,y,z = xyz
        gts.write("   %03i   %+15.8E  %+15.8E  %+15.8E\n"%(atonum,x,y,z))
    gts.write("end_cc\n\n")


    # Write basic data
    rotsigma = "%02i   "%rotsigma
    while len(pgroup) < 5 : pgroup = pgroup+" " 
    gts.write("# Charge, multiplicity, energy [hartree], point group, rotational symmetry number\n")
    gts.write("start_basic\n")
    gts.write("   charge        %i\n"%charge)
    gts.write("   multiplicity  %i\n"%multiplicity)
    gts.write("   energy       %+15.8E # Total energy in hartree\n"%energy)
    gts.write("   pointgroup    %5s          # Point group from automatic assignation\n"%pgroup)
    gts.write("   rotsigma      %5s          # Rotational sigma from automatic assignation\n"%rotsigma)
    gts.write("end_basic\n\n")


    # Write cartesian gradiend
    if g_list is not None and g_list != []:
       gts.write("# Non-scaled cartesian gradient [hartree/bohr]\n")
       gts.write("start_grad\n")
       for gx,gy,gz in g_list:
           gts.write("    %+15.8E  %+15.8E  %+15.8E\n"%(gx,gy,gz))
       gts.write("end_grad\n\n")


    # Write force constant matrix (i.e. hessian matrix)
    if F_list is not None and F_list != []:
       gts.write("# Low triangular part of symmetric force constant (hessian) matrix [hartree/bohr^2]\n")
       gts.write("# i.e.: F_11, F_21, F_22, F_13, F_23, F_33...\n")
       gts.write("start_hess\n")
       for idx in range(0,len(F_list),5):
           line = "  ".join(["%+15.8E"%Fij for Fij in F_list[idx:idx+5]])
           gts.write("           %s\n"%line)
       gts.write("end_hess\n\n")

    # Write freqs
    if F_list is None and freqs is not None:
       gts.write("start_freqs\n")
       for idx in range(0,len(freqs),5):
           line = "  ".join(["%7.2f"%(freq) for freq in freqs[idx:idx+5]])
           gts.write("           %s\n"%line)
       gts.write("end_freqs\n\n")

    gts.close()

def xvecformat(xvec,natoms,out="Nx3"):
    '''
    out: "Nx3" or "3Nx1"
    '''
    # Convert to list in case of numpy array
    try:
       xvec = xvec.tolist()
    except:
       pass

    # xvec in Nx3
    if len(xvec) == natoms:
       if out == "Nx3":
          return np.array(xvec,copy=True)
       if out == "3Nx1":
          final_xvec = []
          for x,y,z in xvec: final_xvec += [x,y,z]
          return np.array(final_xvec,copy=True)
    # xvec in 3Nx1
    elif len(xvec) == 3*natoms:
       if out == "Nx3":
          final_xvec = []
          for idx in range(natoms):
              x,y,z = xvec[3*idx:3*idx+3]
              final_xvec.append( [x,y,z] )
          return np.array(final_xvec,copy=True)
       if out == "3Nx1":
          return np.array(xvec,copy=True)

    sys.exit("Problems in xvecformat function")

def loc1(lista,str1,p1=0):
   pos1 = lista.find(str1,p1)
   return pos1

def parte_str(lista,pos1,pos2):
   trozo = lista[pos1:pos2]
   return trozo

def list_of_lists_to_list(list_of_lists):
    one_list = [item for sublist in list_of_lists for item in sublist] 
    return one_list

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
      #print 'Z-matrix somewhere...'
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
   #print 'Number of atoms = ',natoms
   #print 'Cartesian coordinates'
   #for i in range(natoms):
   #   print atomic_number[i],cartesian[i][0],cartesian[i][1],cartesian[i][2]


   # Merges the list into one element
   length_lista=len(lista_archive)
   lista_archive[0:length_lista]=[''.join(lista_archive[0:length_lista])]
   length_lista=len(lista_archive)
   str_lista=lista_archive[0]
   #print str_lista

   # Command line   
   pcom1=loc1(str_lista,hag,0)
   pcom2=loc1(str_lista,bar2,pcom1)
   commands=parte_str(str_lista,pcom1,pcom2)
   #print 'Command line'
   #print commands

   # Charge and multiplicity
   pchm1=loc1(str_lista,bar2,pcom2+len(bar2))
   pchm2=loc1(str_lista,bar,pchm1+len(bar2))
   t2=parte_str(str_lista,pchm1+len(bar2),pchm2)
   chmul=map(int, t2.split(',') )
   charge=chmul[0]
   mult=chmul[1]
   #print 'Charge = ',charge
   #print 'Multiplicity = ',mult

   # zmatrix or cartesian coordinates
   pgm1=pchm2+1
   pgm2=loc1(str_lista,ver,0)
   #print pgm1,pgm2,str_lista[pgm1],str_lista[pgm2]
   str_geom=parte_str(str_lista,pgm1,pgm2)
   geom=map( str, str_geom.split('\\') )
   if len(geom[0]) <= 4:
   #print str_zmat
      Lzmat=True
      #print 'There is a Z-matrix in the archive file'
      list_zmat=geom
   else:
      Lxyz=True
      xyz_tmp=[]
      tmp_list=[]
      xyz=[]
      ncol_xyz=4
      #print 'There are xyz coordinates in the archive file'
      for i in range(natoms):
	 xyz_tmp.extend(geom[i].split(','))
      for i in range(natoms):
	 atomic_symbol.append(xyz_tmp[ncol_xyz*i])
	 tmp_list=[float(xyz_tmp[i*ncol_xyz+1]),float(xyz_tmp[i*ncol_xyz+2]),float(xyz_tmp[i*ncol_xyz+3])]
	 xyz.append(tmp_list) 
      #for i in range(natoms):
      #	 print atomic_symbol[i],xyz[i][0],xyz[i][1],xyz[i][2]
   if Lxyz:
      ccgeom=xyz
   else:
      ccgeom=cartesian 


   # Energies & other info
   list_level=[]
   list_energy=[]
   #kk=loc1(str_lista,str_en1)
   #if kk < 0: str_first = pgm2
   #else: str_first = str_en1
   pener1x=loc1(str_lista,ver)
   pener1y=loc1(str_lista,str_en1)
   pener1 = max(pener1x,pener1y)
   pener2=loc1(str_lista,bar,pener1)
   pener3=loc1(str_lista,str_en2,pener1)
   t3=parte_str(str_lista,pener2+1,pener3-1)
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
   print 'Energies in Hartree:'
   for i in range(num_energies):
      print 'Energy (',list_level[i],')= ',list_energy[i]

   #Hessian
   if Lhess:
      phes1=loc1(str_lista,str_imag,0)
      phes2=loc1(str_lista,bar2,phes1)
      phes3=loc1(str_lista,bar2,phes2+len(bar2))
      #print str_lista[phes1],str_lista[phes2],str_lista[phes3]
      num_imag=int(parte_str(str_lista,phes1+len(str_imag),phes2))
      #print 'Number of imaginary frequencies',num_imag
      str_hess=parte_str(str_lista,phes2+len(bar2),phes3)
      #print str_hess
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
      #print 'No hessian in this file'
      pass

   #Gradient
   if Lhess > 0:
      pgrd1=phes3+len(bar2)
      pgrd2=loc1(str_lista,fin,0)
      str_grd=parte_str(str_lista,pgrd1,pgrd2-3)
      grad=map( float, str_grd.split(',') )
      #print grad
   else:
      #print 'No gradient in this file'
      pass
   
   if not Lzmat: list_zmat=[0] 
   if not Lhess:
      hess_data=[0]
      hess_full=[0]
      grad=[0]

   return num_imag,charge,mult,commands,list_level,list_energy,atomic_number,ccgeom,\
          list_zmat,hess_data,hess_full,grad

def f_mstor_gen(list_zmat,natoms):
    ##ZMAT for dihedral automatic search
    ##writes the results in f_tor
    lz_mod = []
    print list_zmat
    target_idx = natoms + 2
    lz_mod = list_zmat[:target_idx] 
    print lz_mod
    ntor_scan = 0
    ntor_ch3 = 0
    pos_tor = []
    pos_ch3 = []
    lsplt_tor = []
    lsplt_ch3 = [] 
    list_tor = []
    list_ch3 = []
    list_sigtor = []
    idx = 0
    print lz_mod
    for str1 in lz_mod:
        if 'tor' in str1:
            ntor_scan   += 1
            pos_tor.append(idx)
        if 'CH3' in str1:
            ntor_ch3    += 1
            pos_ch3.append(idx)
        idx += 1
    print pos_tor,pos_ch3
    for t in range(ntor_scan):
        str_tor = lz_mod[pos_tor[t]]
        lsplt_tor = str_tor.split(',')   
        tup_tor=(pos_tor[t]+1,int(lsplt_tor[1]),int(lsplt_tor[3]),int(lsplt_tor[5]))
        list_tor.append(tup_tor)
        list_sigtor.append(1)
    if ntor_ch3 > 0:
        for t in range(ntor_ch3):
            str_ch3 = lz_mod[pos_ch3[t]]
            lsplt_ch3 = str_ch3.split(',')   
            tup_ch3=(pos_ch3[t]+1,int(lsplt_ch3[1]),int(lsplt_ch3[3]),int(lsplt_ch3[5]))
            list_ch3.append(tup_ch3)
            list_sigtor.append(3)
    l_tor_tot = list_tor + list_ch3
    n_tor_tot = ntor_scan + ntor_ch3
    ## Print
    #f = open(f_tor,'w+')
    #f.write(" %3i\n"%n_tor_tot)
    #for i in range(n_tor_tot):
    #   f.write(('{:3d}'*4+"\n").format(*l_tor_tot[i]))
    #f.write('\n')
    #f.write(('{:3d}'*n_tor_tot+"\n").format(*list_sigtor))
    #f.write('\n')
    #f.close()
    return n_tor_tot,l_tor_tot,list_sigtor

def main(filename,print_option=None):

   # Atomic masses
   mass_dict = {1:1.007825e+00,2:4.0026e+00,3:7.01600e+00,4:9.01218e+00,5:11.00931e+00,
        6:12.0e+00,7:14.00307e+00,8:15.99491e+00,9:18.99840e+00,10:19.99244e+00,
        11:22.9898e+00,12:23.98504e+00,13:26.98153e+00,14:27.97693e+00,15:30.97376e+00,
        16:31.97207e+00,17:34.96885e+00,18:39.948e+00,19:38.96371e+00,20:39.96259e+00,
        21:44.95592e+00,22:47.90e+00,23:50.9440e+00,24:51.9405e+00,25:54.9381e+00,
        26:55.9349e+00,27:58.9332e+00,28:57.9353e+00,29:62.9298e+00,30:63.9291e+00,
        31:68.9257e+00,32:73.9219e+00,33:74.9216e+00,34:79.9165e+00,35:78.9183e+00,
        36:83.9115e+00,37:84.9117e+00,38:87.9056e+00,39:89.9054e+00,40:89.9043e+00,
        41:92.9060e+00,42:97.9055e+00,43:97.0e+00,44:101.9037e+00,45:102.9048e+00,
        46:105.9032e+00,47:106.9041e+00,48:113.9036e+00,49:114.9041e+00,50:119.9022e+00,
        51:120.9038e+00,52:129.9067e+00,53:126.9044e+00,54:131.9042e+00,55:132.9054e+00,
        56:137.9052e+00,57:138.9063e+00,58:139.9054e+00,59:140.9076e+00,60:141.9077e+00,
        61:144.9127e+00,62:151.9197e+00,63:152.9212e+00,64:157.9241e+00,65:158.9253e+00,
        66:163.9292e+00,67:164.9303e+00,68:165.9303e+00,69:168.9342e+00,70:173.9389e+00,
        71:174.9408e+00,72:179.9465e+00,73:180.9480e+00,74:183.9509e+00,75:186.9557e+00,
        76:191.9615e+00,77:192.9629e+00,78:194.9648e+00,79:196.9665e+00,80:201.9706e+00,
        81:204.9744e+00,82:207.9766e+00,83:208.9804e+00,84:208.9824e+00,85:209.9871e+00,
        86:222.0176e+00,87:223.0197e+00,88:226.0254e+00,89:227.0278e+00,90:232.0381e+00,
        91:231.0359e+00,92:238.0508e+00,93:237.0482e+00,94:244.0642e+00,95:243.0614e+00,
        96:247.0703e+00,97:247.0703e+00,98:251.0796e+00,99:252.0829e+00,100:257.0751e+00,
        101:258.0986e+00,102:259.1009e+00,103:260.1053e+00}

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
   print 'There are ',number_archives,' archive files in this Gaussian file'
   if number_archives == 1 and start_archives[0] == 0:
      print 'This Gaussian file does not have an archive part'
      exit()

   str_geom1s=[]
   str_geom2s=[]
   str_nosyms=[]
   str_zmats=[]
   for i in range(number_archives):
      print '---------------------------'
      print 'Archive file number : ',i
      print '---------------------------'
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


      # ------------ For MSTor --------------
      # With the first electronic structure output file:
      #   i)  Generates the torsions and store them in a list
      if print_option == 20:
          if i == 0 and list_zmat != []:
                n_tor_tot,l_tor_tot,l_sigtor = f_mstor_gen(list_zmat,natoms) 
          print natoms,mult,n_tor_tot,l_tor_tot,l_sigtor 
          return natoms,mult,n_tor_tot,l_tor_tot,l_sigtor


      # ------------ END For MSTor --------------
         
      #print 'Commands of Gaussian input'
      #print commands
      #print 'List of levels'
      #print list_level
      #print 'List of energies'
      #print list_energy
      #print 'Charge = ',charge
      #print 'Multiplicity = ',mult
      #print 'Atomic numbers'
      #print atomic_number
      #print 'Cartesian coordinates geometry'
      #print ccgeom
      #print 'Updated Z-matrix'
      #for i in range(len(list_zmat)):
      #   print list_zmat[i]
      #print 'Number of imaginary frequencies'
      #print num_imag
      #print 'Hessian'
      #print hess_data
      #for idx in range(0,len(hess_data),6):
      #   line="".join([" %+11.8f"%Fij for Fij in hess_data[idx:idx+6]])
      #   print line
      #print 'Gradient'
      #print grad

   # atomic mass
   atomic_mass = []
   for i in range(natoms):
      atomic_mass.append(mass_dict[atomic_number[i]])

   if idx_last_zmat > -1 and idx_last_zmat < number_archives-1 and print_option == 11:
      print 'WARNING: The printed Z-matrix is not from the last archive file'
      list_zmat = list_zmat_backup[:]

   # point group symmetry and symmetry number of rotation
   pgroup,rotsigma = get_pgs(atomic_number,atomic_mass,ccgeom)

   if print_option == None:
      return num_imag
   elif print_option < 50:
      file_prt="prt_archive"
      fprt=open(file_prt,'w')
      ccgeom1 = list_of_lists_to_list(ccgeom)
      if print_option == 0:
         fprt.write("%s"%str(num_imag))
      if print_option == 1: # write the geometry and the force constants
         fprt.write(' %s,%s\n'%(str(charge),str(mult)))
         lines_list_zmat = len(list_zmat)
         if lines_list_zmat > 1:
            for i in range(lines_list_zmat-1):
               if i == lines_list_zmat-2 and list_zmat[i] == '':
                  pass
               else: 
                  fprt.write('%s\n'%list_zmat[i])
            fprt.write('\n')
         else:
            for idx in range(0,len(ccgeom1),3):
               x,y,z = ccgeom1[idx+0:idx+3]
               i=idx/3
               atonum = atomic_number[i]
               fprt.write("   %3i   %+15.8E  %+15.8E  %+15.8E\n"%(atonum,x,y,z))
            fprt.write('\n')
         fprt.write(' %+11.8f\n'%0.)
         for idx in range(0,len(grad),6):
            line="".join([" %+11.8f"%gij for gij in grad[idx:idx+6]])
            fprt.write(line+"\n")
         for idx in range(0,len(hess_data),6):
            line="".join([" %+11.8f"%Fij for Fij in hess_data[idx:idx+6]])
            fprt.write(line+"\n")
      fprt.close() 

      if print_option == 2:  # write gts file
         energy = list_energy[-1]
         ccgeom1 = list_of_lists_to_list(ccgeom)
         # from Angstrom to Bohr
         for i in range(3*natoms):
            ccgeom1[i] = ccgeom1[i]/0.52917721067
         write_gtsfile(ccgeom1 , atomic_number , charge, mult, energy, \
                   pgroup, rotsigma , grad, hess_data, file_prt) 

      if print_option == 11: # write the geometry and the force constants
         file_prtg="prt_geom"
         fprtg=open(file_prtg,'w')
         #fprtg.write(' %s,%s\n'%(str(charge),str(mult)))
         lines_list_zmat = len(list_zmat)
         if lines_list_zmat > 1:
            for i in range(lines_list_zmat):
               if i == lines_list_zmat-2 and list_zmat[i] == '':
                  pass
               elif i == lines_list_zmat-1 and list_zmat[i] == '':
                  pass
               elif i == lines_list_zmat-1 and list_zmat[i] == 'Test job not archived.':
                  pass
               elif list_zmat[i].lower() == 'variables:':
                  fprtg.write('\n')
               else:    
                  fprtg.write('%s\n'%list_zmat[i].replace(',',' '))
         else:
            for idx in range(0,len(ccgeom1),3):
               x,y,z = ccgeom1[idx+0:idx+3]
               i=idx/3
               atonum = atomic_number[i]
               fprtg.write("   %3i   %+15.8E  %+15.8E  %+15.8E\n"%(atonum,x,y,z))
         fprtg.close()
         file_prth="prt_hess"
         fprth=open(file_prth,'w')
         fprth.write(' %+11.8f\n'%0.)
         for idx in range(0,len(grad),6):
            line="".join([" %+11.8f"%gij for gij in grad[idx:idx+6]])
            fprth.write(line+"\n")
         for idx in range(0,len(hess_data),6):
            line="".join([" %+11.8f"%Fij for Fij in hess_data[idx:idx+6]])
            fprth.write(line+"\n")
         fprth.close()
   else:
      if print_option == 99:  # 
         print pgroup,rotsigma 
         return pgroup,rotsigma
      if print_option == 101:
         energy = list_energy[-1]
         return energy
      return 

   
# This part is executed only if archive.py is the main program
if __name__ == "__main__":

    print main("Z.log",print_option=101)
    exit()
    from sys import argv
    args = sys.argv[1:]

    filename   = args[0]

    # options
    if "--prt"    in args:
       ii  = args.index("--prt")
       print_option = int(args[ii+1])
    else:
       print_option = None

    if print_option == None:
       num_imag=main(filename)
    else:
       main(filename,print_option)
  
