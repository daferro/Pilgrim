# _Pilgrim_

## About _Pilgrim_

    Name of the Program: Pilgrim
    Program Version : 1.0  
    Program Version Date: April 3, 2019  
    Manual  Version Date: April 3, 2019

_Pilgrim_ is an user-friendly program written in Python 2.7.  
It was designed to use direct-dynamics to calculate thermal rate constants of 
chemical reactions and to simulate chemical kinetics mechanisms.

For reaction processes with many elementary steps, each of these steps can be 
calculated using conventional transition state theory (TST) or variational transition 
state theory (VTST). In this version, _Pilgrim_ can calculate thermal rate constants 
with the canonical version of the variational transition state theory (CVT), 
which requires the calculation of the minimum energy path (MEP) associated to each elementary step.
Moreover, multi-dimensional quantum effects can be incorporated through the 
small-curvature tunneling approximation (SCT). 
The above methodologies are available for reactions involving a single 
structure and for reactions involving flexible molecules with multiple conformations. 
Specifically, for systems with many conformers the program can evaluate each of 
the elementary reactions by multi-path canonical variational transition 
state theory (MP-CVT) or multi-structural VTST (MS-VTST). 
Torsional anharmonicity can be also incorporated through the MSTor and Q2DTor programs.

_Pilgrim_ also performs dual-level calculations automatically. 
First, low-level calculations are carried out for the reaction of interest and second, 
single-point energy calculations of the reactants, transition state, points along the MEP 
and products are performed at a higher level. 
The low-level calculations are corrected with the high-level single point energies 
using the interpolated single-point energies (ISPE) algorithm. 

Once all the rate constants of the chemical processes of interest are known, 
by means of their calculation using _Pilgrim_ or by using an analytical expression, 
it is possible to simulate the whole process using kinetic Monte Carlo (KMC). 
This algorithm allows performing a kinetics simulation and monitoring the evolution 
of each chemical species with time, as well as providing its chemical yield. 


## How to cite

D. Ferro-Costas, D. G. Truhlar, A. Fernández-Ramos, Pilgrim - version 1.0
(University of Minneapolis, Minnesota, MN, and Universidade de Santiago
de Compostela, Spain, 2019). https://github.com/daferro/Pilgrim


## Licensing and Distribution 

_Pilgrim version 1.0_

Copyright (c) 2019, David Ferro Costas (david.ferro@usc.es) and Antonio Fernandez Ramos (qf.ramos@usc.es)

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


## Description of files

 Contents of the folders distributed in this version:
  - **src/**       : _Pilgrim_ source files
  - **docs/**      : Manual of _Pilgrim_
  - **tests/**     : All the files related to the tests set
        

## Installation

_Pilgrim_ is a program written in Python 2. Consequently, it does not need any kind 
of compilation, as it would be the case with C or Fortran programs.
The user should install Python 2 in order to use _Pilgrim_ (version 2.7 is recommended), 
as well as the following Python libraries:
   - cmath
   - fcntl
   - glob
   - math
   - matplotlib
   - multiprocessing
   - numpy
   - os
   - random
   - scipy
   - sys
   - time

WARNING: __do not__ use Python 3 to execute _Pilgrim_.


## Setting up the program

Before using _Pilgrim_, the user has to define the path to the executable(s) of the 
software for the electronic structure calculation (ESSO).
At the moment _Pilgrim_ supports __Gaussian__ and __Orca__ as ESSO.

In order to interact with the ESSO, _Pilgrim_ needs to know the location of some executable files. 
Such information is obtained from the following environment variables, which have to be 
defined and exported by the user in their __.bashrc__ file:

For __Gaussian__ users, the environment variables are:

  - _GauExe_, the path to the Gaussian executable and
  - _GauFchk_, path to the formchk tool.

Example of paths are:

  - ```export Gauexe="/home/programs/G09_64D/g09/g09"``` 
  - ```export GauFchk="/home/programs/G09_64D/g09/formchk"```  

Similarly, for __Orca__ users, this file has to contain the environment variable:

  - _OrcaExe_, the path to the Orca executable.

Example:

  - ```export OrcaExe="/home/programs/orca_4_0_1_2/orca"```  

where, again, the path to the Orca executable is between quotation marks. 


## Execution

You can run _Pilgrim_ by invoking the Python interpreter manually as follows:

```python2.7 pilgrim.py```

If you prefer to avoid invoking the Python interpreter, you have to follow these two simple steps:

Add as the first line in the pilgrim.py file the following:

```#!PATH_FOR_PYTHON python```

where PATH_FOR_PYTHON indicates the location of the Python interpreter.

Example:

```#!/usr/bin/env python```

In this example Python is located in `/usr/bin/env`.

Make the main program pilgrim.py executable:

```chmod u+x pilgrim.py```

This allows you to run _Pilgrim_ just using:

```pilgrim.py```

Before run _Pilgrim_, we recommend to read the help menu. It can be displayed either by typing

```pilgrim.py --help```

or

```pilgrim.py -h```

## Tests set

Directory __tests/__ contains two subdirectories:

  - __tests\_to\_run/__:   
    Includes the electronic structure files (ESFILs) with the 
    stationary points needed to run each of the working examples (WEs) 
    using _Pilgrim_. See the manual for details.
     
  - __tests\_out/__:   
    Includes the results of running the WEs. They can also be used
    for comparison.


                                                            

