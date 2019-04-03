#!/usr/bin/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.checkmods  |
| Last Update:  2019/01/15 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*

The function in this module is in
charge of asserting that the Python
libraries imported by Pilgrim are
actually installed
'''

#===============================================================#
def checkmods():
    failed  = []
    modules = ["cmath","fcntl","gc","glob","math","matplotlib",\
               "multiprocessing","numpy","os","random","scipy",\
               "sys","time"]
    for module in modules:
        try:
          # because we want to import using a variable, do it this way
          module_obj = __import__(module)
          # create a global object containing our module
          globals()[module] = module_obj
        except ImportError:
          failed.append(module)
    if failed != []:
       print "ERROR! Missing Python module(s):"
       for module in failed: print "    * %s"%module
       print "Install them before executing Pilgrim!"
       exit(1)
#===============================================================#

