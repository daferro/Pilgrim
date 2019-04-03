# level: RHF STO-3G
# Atomic number and non-scaled cartesian coordinates [bohr]
start_cc
   001   +0.00000000E+00  +0.00000000E+00  +2.62322829E-02
   001   +0.00000000E+00  +0.00000000E+00  +1.37216506E+00
end_cc

# Charge, multiplicity, energy [hartree],
# point group and rotational symmetry number
start_basic
   charge        0
   multiplicity  1
   energy       -1.11750589      # Total energy in hartree
   pointgroup    Dinfv           # Point group
   rotsigma      2               # Rotational sigma
end_basic

# Non-scaled cartesian gradient [hartree/bohr]
start_grad
   +0.00000000E+00  +0.00000000E+00  -7.68200073E-06
   +0.00000000E+00  +0.00000000E+00  +7.68200073E-06
end_grad

# Low triangular part of force constant (hessian) matrix [hartree/bohr^2]
# i.e.: F_11, F_21, F_22, F_13, F_23, F_33...
start_hess
   +5.70756646E-06  +0.00000000E+00  +5.70756646E-06  +0.00000000E+00  +0.00000000E+00
   +5.72909758E-01  -5.70756646E-06  +0.00000000E+00  +0.00000000E+00  +5.70756647E-06
   +0.00000000E+00  -5.70756647E-06  +0.00000000E+00  +0.00000000E+00  +5.70756647E-06
   +0.00000000E+00  +0.00000000E+00  -5.72909758E-01  +0.00000000E+00  +0.00000000E+00
   +5.72909758E-01
end_hess

