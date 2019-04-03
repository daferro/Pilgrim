# |sin(x)|   < EPS_SCX --> sin(x) =  0
# |cos(x)-1| < EPS_SCX --> cos(x) =  1
# |cos(x)+1| < EPS_SCX --> cos(x) = -1
EPS_SCX = 1e-7

# small angle
EPS_SMALLANGLE = 1e-3

# ds STEP FOR THE CALCULATION OF RETURN POINTS
DS_RPT = 1e-3

# comparision of geometries (bohr)
EPS_GEOM = 1e-5

# zero moment of inertia
EPS_INERTIA = 1e-7

# the norm of a vector is ZERO
EPS_NORM = 1e-7

# step for numerical calculation of hessian matrix in analitic surfaces
EPS_HESSDX = 1e-4

# step for cubic first step in MEP calculation 
STEP_CUBIC = 1e-4

# comparision of float numbers
EPS_FLOAT = 1e-8

# comparision of float numbers (strict)
EPS_FSTRICT = 1e-10

# comparision of temperatures
EPS_TEMP = 1e-8

# zero por single-value decomposition
EPS_SVD = 1e-9

# zero por generalized inverse
EPS_GIV = 1e-10

# Ignore ic-freqs smaller than EPS_ICF cm^-1
EPS_ICF = 0.3

# If cc-freq and ic-freq differ less than EPS_CCIC cm^-1, they can be considered to be equal
EPS_CCIC = 0.5

# comparison of s values of MEP (in bohr)
EPS_MEPS = 1e-7

# comparison of energies of MEP (in hartree)
EPS_MEPE = 1e-8

# to identify s values in DLEVEL
EPS_DLEVELS = 1e-9

# zero for the comparison of masses in amu
EPS_AMU = 0.001

