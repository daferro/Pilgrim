import numpy as np
import random
import time
import physcons as pc

from Exceptions import NoReacMol
from criteria import EPS_FSTRICT

#---------------------
def read_kmc_inp(ifile):
    '''
    Input file, AFR version
    '''
    # Some initializations
    rk=[]
    s_init1=[]
    s_init2=[]
    s_fin1=[] 
    s_fin2=[] 
    xm0=[]
    mol_label=[]
    reac_ex=0
    xm0cons=0.

    inputkmc=open(ifile,'r') 
    titulo=inputkmc.readline()

    # Number of species (nsp)  and initial number of molecules
    # Excess reactant has the string 'excess' at the end (constant concentration)
    nsp=int(inputkmc.readline())
    iexces=0
    for idum in range(nsp):
       leelinea=inputkmc.readline()
       if 'excess' in leelinea:
          kk1,x1,_label,_labelx=leelinea.split()
          reac_ex +=1
          iexces=idum+1
       else:
          kk1,x1,_label=leelinea.split()
       xm0.append(float(x1))
       mol_label.append(str(_label))
    if reac_ex > 1:
       print 'Only one species in excess can be considered'
       exit()
    if reac_ex == 1:
        xm0cons=xm0[iexces-1]
    
    # rate constants at a given temperature
    # species for the initial and final states
    # deals with unimolecular and bimolecular reactions  
    nrx=int(inputkmc.readline())
    for idum in range(nrx):
       leelinea=inputkmc.readline()
       rklee,s_i_lee1,s_i_lee2,s_f_lee1,s_f_lee2=leelinea.split()
       rk.append(float(rklee))
       s_init1.append(int(s_i_lee1))
       s_init2.append(int(s_i_lee2))
       s_fin1.append(int(s_f_lee1))
       s_fin2.append(int(s_f_lee2))
    # Plotting
    splot=float(inputkmc.readline())
    # volume (reads the volume if there are bimolecular rate constants) 
    s_bim = sum([s_init1[i]*s_init2[i] for i in range(nrx)])
    if s_bim > 0:
       volume = float(inputkmc.readline())
    else:
       volume = 1. 
    # closes file
    inputkmc.close()

    # Convert data to good format
    print "s_init1:", s_init1
    print "s_init2:", s_init2
    print "s_fin1:", s_fin1
    print "s_fin2:", s_fin2
    print "mol_label:", mol_label
    print "xm0cons:", xm0cons
    print "xm0:", xm0
    print "reac_ex:", reac_ex
    print "rk:", rk
    print "splot:", splot
    print "volume:", volume
    print "s_bim:", s_bim
    print "titulo:", titulo
    print "iexces:", iexces

    processes = {}
    for idx in range(nrx):
        R1 = s_init1[idx]
        R2 = s_init2[idx]
        P1 = s_fin1[idx]
        P2 = s_fin2[idx]
        if R1 != 0: R1 = mol_label[R1-1]
        else      : R1 = None
        if R2 != 0: R2 = mol_label[R2-1]
        else      : R2 = None
        if P1 != 0: P1 = mol_label[P1-1]
        else      : P1 = None
        if P2 != 0: P2 = mol_label[P2-1]
        else      : P2 = None
        if None in [R1,R2]: k_secs = rk[idx]
        else              : k_secs = rk[idx] / volume
        processes[idx+1] = (rk[idx],k_secs,R1,R2,P1,P2)
    
        
    diconcs = {}
    for idx in range(len(xm0)):
        label = mol_label[idx]
        diconcs[label] = xm0[idx]

    if iexces != 0:
       excess_specie = mol_label[iexces-1]
    else:
       excess_specie = None

    volume = volume
    nstp   = splot

    # Return data
    return diconcs, processes, excess_specie, volume, nstp
#---------------------


#=====================================================#
def calculate_propensities(dxvec, processes):
    '''
    Calculates the propensitie of each reaction
    '''
    propensities = []
    total_propensity = 0.0
    for Rs,Ps,ks in processes:
        prop = ks
        for Ri in Rs: prop *= dxvec[Ri]
        propensities.append( prop )
        total_propensity += prop
    return propensities, total_propensity
#-----------------------------------------------------#
def generate_random(eps=EPS_FSTRICT):
    '''
    Generates a random number, excluding those smaller than eps
    '''
    while True:
      randx=random.random() 
      if randx >= eps: return randx
#=====================================================#




#=====================================================#
def kmc(ipops, processes, excess_species=None, volume=1.0/pc.ML, nstpdata=1000):
    '''
    ipops   : dictionary with initial populations (only needed those != 0.0)
    processes: a list with the reactants, products and the rate constant
               processes[idx] = (Rs,Ps,k)

    INPUT UNITS: atomic units
                 concentrations: molecules/bohr**3
                 volume        : bohr**3
                 rate constants: in au
    '''

    #----------------------------#
    # set initial concentrations #
    # and reactant molecules     #
    #----------------------------#
    reactants = []
    for Rs,Ps,k in processes:
        for species in Rs+Ps:
            if species not in ipops.keys(): ipops[species] = 0.0
            bool1 = species in excess_species
            bool2 = species in reactants
            if (not bool1) and (not bool2): reactants.append(species)

    #------------------------#
    # rate constants to s^-1 #
    #------------------------#
    for idx,(Rs,Ps,k) in enumerate(processes):
        nR = len(Rs)
        k /= volume**(nR-1)
        processes[idx] = (Rs,Ps,k)

    #-----------------------------------------#
    # Get dict of xvec and limiting molecules #
    #-----------------------------------------#
    dxvec = ipops.copy()
    try   : N0 = min([pop for pop in dxvec.values() if pop != 0.0])
    except: N0 = 0.0
    if N0 == 0.0: raise NoReacMol(Exception)
        
    # --------------------------#
    # START KINETIC MONTE-CARLO #
    # --------------------------#
    # initialize variables
    tau, tx = 0.0, 0.0
    jcount  = 0

    # data
    xvalues = [tx]
    yvalues = {key:[val] for key,val in dxvec.items()}

    Nj = N0
    xi = np.array([ dxvec[species] for species in sorted(dxvec.keys())])
    while Nj > 0.0:
       # compare each 1000 steps
       if jcount % 1000 == 0 and jcount > 0:
          xj = np.array([ dxvec[species] for species in sorted(dxvec.keys())])
          diff = np.linalg.norm(xj-xi)/1000
          if diff < 1e-4: break
          xi = xj
       # calculate propensities
       propensities, tot_propensity = calculate_propensities(dxvec,processes)
       #if tot_propensity < 1.e-15: break
       # generate two random numbers 
       r1=generate_random()    # 1) Tau (time)
       r2=generate_random()    # 2) Changes in population
       # Select process (stacking)
       value = tot_propensity * r2
       sum_props = 0.0
       for target,prop in enumerate(propensities):
           sum_props += prop
           if sum_props >= value: break
       # Modify populations (except those in excess)
       Rs,Ps,ks = processes[target]
       for Ri in Rs:
           if Ri not in excess_species: dxvec[Ri] -= 1.0
       for Pi in Ps:
           if Pi not in excess_species: dxvec[Pi] += 1.0
       # Time step
       if tot_propensity == 0.0: break
       tau = np.log(1./r1)/tot_propensity
       tx += tau
       # Keep data
       jcount += 1
       if jcount%nstpdata==0:
         xvalues.append(tx)
         for specie in dxvec.keys():
             yvalues[specie].append( dxvec[specie] )
    
       # Calculate current number of reactant molecules
       Nj = sum([dxvec[specie] for specie in reactants])

    return xvalues, yvalues
#=====================================================#


if __name__ == '__main__':
   from sys import argv
   script, fkmcin = argv
   diconcs, processes, excess_specie, volume, nstp = read_kmc_inp(fkmcin)
   kmc(diconcs, processes, excess_specie, volume, nstp)

