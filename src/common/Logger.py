#!/usr/bin/python2.7
'''
*----------------------------------*
| Package    :  common             |
| Module     :  Logger             |
| Last Update:  2018/09/04 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the Logger class
'''

import sys

class Logger(object):
    '''
    Class used to save in a file
    what is printed in terminal
    Requirements:
      sys library
    Use:
      sys.stdout = Logger(f_out,tprint)
    '''

    def __init__(self,output=None,mode="w",bool_print=True):
        # terminal
        self.terminal   = sys.__stdout__
        self.bool_print = bool_print

        # file
        self.bool_write = False
        self.file       = output
        if output is not None:
           self.log = open(self.file, mode)
           self.bool_write = True
        else:
           self.bool_write = False

    def write(self, message):
        if self.bool_print: self.terminal.write(message)
        if self.bool_write: self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

