#!/usr/bin/python2.7
'''
*-----------------------------------*
| Last Update:  2019/04/03 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*
'''

PROGNAME  = '''
 -----------------------------------------------------------
  __________.___.____     __________________.___   _____    
  \______   \   |    |   /  _____/\______   \   | /     \   
   |     ___/   |    |  /   \  ___ |       _/   |/  \ /  \  
   |    |   |   |    |__\    \_\  \|    |   \   /    Y    \ 
   |____|   |___|_______ \______  /|____|_  /___\____|__  / 
                        \/      \/        \/            \/  
 -----------------------------------------------------------\n'''

VERSION = "1.0 (2019-04-03)"
AUTHORINFO = '''
 -----------------------------------------------------------
  Program version: Pilgrim %s
 -----------------------------------------------------------
  A thermal rate constant calculator and kinetics Monte
  Carlo Simulator.
                                                         
  Main authors:                                              
       * Ferro-Costas, David       (1)                       
       * Fernandez-Ramos, Antonio  (1)                       
                                                          
  In collaboration with:                                     
       * Truhlar, Donald G.        (2)                       
                                                          
  (1) Centro Singular de Investigacion en Quimica Bioloxica
      e Materiais Moleculares (CIQUS), Universidade de
      Santiago de Compostela, Galicia, Spain 

  (2) Department of Chemistry and Supercomputer Institute,
      University of Minnesota, Minneapolis, Minnesota        
 -----------------------------------------------------------\n'''%VERSION


#==========================================================#
# FIRSTLY: check modules!!!                                #
#----------------------------------------------------------#
import mpilgrim.checkmods   as     checkmods               #
checkmods.checkmods()                                      #
#----------------------------------------------------------#
import os                                                  #
import sys                                                 #
import time                                                #
#----------------------------------------------------------#
import common.Exceptions       as     Exc                  #
from   common.fncs             import classify_args        #
from   common.fncs             import print_string         #
from   common.fncs             import time2human           #
#----------------------------------------------------------#
import mpilgrim.names       as PN                          #
import mpilgrim.rdwr        as RW                          #
#----------------------------------------------------------#
import mpilgrim.opt_gather     as gather                   #
import mpilgrim.opt_input      as inpmenu                  #
import mpilgrim.opt_ics        as ics                      #
import mpilgrim.opt_pfn        as pfn                      #
import mpilgrim.opt_path       as path                     #
import mpilgrim.opt_rcons      as rcons                    #
import mpilgrim.opt_kmc        as kmc                      #
import mpilgrim.opt_fit        as fit                      #
import mpilgrim.opt_plot       as plot                     #
import mpilgrim.opt_hlcalc     as hlcalc                   #
from   mpilgrim.exceptions     import deal_with_exception  #
from   mpilgrim.diverse        import get_input_data       #
from   mpilgrim.diverse        import dlevel_to_files      #
#----------------------------------------------------------#
from   mpilgrim.helps          import HELP_main            #
from   mpilgrim.helps          import HELP_gather          #
from   mpilgrim.helps          import HELP_input           #
from   mpilgrim.helps          import HELP_ics             #
from   mpilgrim.helps          import HELP_pfn             #
from   mpilgrim.helps          import HELP_path            #
from   mpilgrim.helps          import HELP_hlcalc          #
from   mpilgrim.helps          import HELP_rcons           #
from   mpilgrim.helps          import HELP_fit             #
from   mpilgrim.helps          import HELP_kmc             #
from   mpilgrim.helps          import HELP_plot            #
#==========================================================#


#==========================================================#
#                    Checking arguments                    #
#==========================================================#
def check_ics(dargs):
    key = "ics"
    if key in dargs.keys():
        if   len(dargs[key]) == 0:
             dargs[key].append("2")
             dargs[key].append("*")
        elif len(dargs[key]) == 1:
             dargs[key].append("*")
        elif len(dargs[key]) == 2:
             pass
        else:
            print "Maximum number of variables for --ics is two!!"
            exit()
    return dargs
#----------------------------------------------------------#
def check_sft(dargs):
    key = "software"
    if key in dargs.keys():
        try   : return dargs[key][0]
        except: return None
    else      : return "gaussian"
#----------------------------------------------------------#
def check_dlevel(dargs):
    key = "dlevel"
    if key in dargs.keys(): return True
    else                  : return False
#==========================================================#



#==========================================================#

def main():
    date = time.strftime("%Y/%m/%d %H:%M")
    t1 = time.time()

    # Read user arguments
    user_args = sys.argv[1:]
    if len(user_args) == 0:
       print PROGNAME
       print AUTHORINFO
       print_string(HELP_main,2)
       return
    dargs = classify_args(user_args)

    # User asked for version
    if   "version"    in dargs.keys():
       print "Current version: %s"%VERSION
       return

    # User asked for help
    elif "help"    in dargs.keys():
       print PROGNAME
       print AUTHORINFO
       if len(dargs["help"]) == 0   : print_string(HELP_main  ,2)
       if "gather"  in dargs["help"]: print_string(HELP_gather,1)
       if "input"   in dargs["help"]: print_string(HELP_input ,1)
       if "ics"     in dargs["help"]: print_string(HELP_ics   ,1)
       if "pfn"     in dargs["help"]: print_string(HELP_pfn   ,1)
       if "path"    in dargs["help"]: print_string(HELP_path  ,1)
       if "hlcalc"  in dargs["help"]: print_string(HELP_hlcalc,1)
       if "rcons"   in dargs["help"]: print_string(HELP_rcons ,1)
       if "fit"     in dargs["help"]: print_string(HELP_fit   ,1)
       if "kmc"     in dargs["help"]: print_string(HELP_kmc   ,1)
       if "plot"    in dargs["help"]: print_string(HELP_plot  ,1)
       return

    # Print logo and current date
    print PROGNAME
    print AUTHORINFO
    print
    print " -----------------------------------------------------------"
    print "  Current date: %s"%date
    print " -----------------------------------------------------------"
    print 

    # check some arguments
    dargs    = check_ics(dargs)
    software = check_sft(dargs)
    dlevel   = check_dlevel(dargs)

    # --software used properly
    if software is None:
       print "  Software option has to be followed by an argument!"
       print
       exit()

    # Act according user argument
    IN_OPTS = False
    for option in "ls,gather,input,ics,pfn,path,hlcalc,rcons,kmc,fit,plot".split(","):
        if option not in dargs.keys(): continue
        IN_OPTS = True
        # read input files and print table
        idata, status, string = get_input_data()
        print_string(string,nbs=3)
        # data files according to case
        datafiles,string = dlevel_to_files(dlevel)
        if option in "pfn,path,hlcalc,rcons,kmc,fit,plot".split(","):
           print_string(string,nbs=3)
        case = (datafiles,dlevel,software)
        # execute pilgrim with option
        print  " ==========================================================="
        print  " ||  EXECUTING PILGRIM WITH --%-6s                      ||"%option
        print  " ==========================================================="
        print
        if   option == "ls"    :  gather.ls_struc(idata[0][0])
        elif option == "gather":  gather.main(idata,status                      )
        elif option == "input" : inpmenu.main(idata,status                      )
        elif option == "ics"   :     ics.main(idata,status,*dargs["ics"   ]     )
        elif option == "pfn"   :     pfn.main(idata,status,case, dargs["pfn"   ])
        elif option == "path"  :    path.main(idata,status,case, dargs["path"  ])
        elif option == "hlcalc":  hlcalc.main(idata,status,case, dargs["hlcalc"])
        elif option == "rcons" :   rcons.main(idata,status,case, dargs["rcons" ])
        elif option == "kmc"   :     kmc.main(idata,status,case                 )
        elif option == "fit"   :     fit.main(idata,status,case, dargs["fit"   ])
        elif option == "plot"  :    plot.main(             case, dargs["plot"  ])
        print
        print " ===========================================================\n"
        break

    # User asked for nothing. Print help!
    if not IN_OPTS: print_string(HELP_main,2); return


    # Print elapsed time
    t2 = time.time()
    timeline = "Total Elapsed Time: %5.1f %5s |"%time2human(t2-t1,"secs")
    print "                            "+timeline
    print "                            "+"-"*len(timeline)
#==========================================================#


#==========================================================#
if __name__ == '__main__':
   try: main()
   except Exception as exception: deal_with_exception(exception)
#==========================================================#


