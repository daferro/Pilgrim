#!/usr/bin/python2.7
'''
*----------------------------------*
| Package    :  common             |
| Module     :  MyCompleter        |
| Last Update:  2018/09/04 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the MyCompleter class
as well as a function to set it
'''

import glob

class MyCompleter(object):
      '''
      Custom completer
      '''
      def __init__(self, options=[], sort=False):
          self.options = options
          if sort: self.options.sort()
      def complete(self, text, state):
          # on first trigger, build possible matches
          if state == 0:
              # cache matches (entries that start with entered text)
              if text:
                  self.matches = [s for s in self.options if s and s.startswith(text)]
              # no text entered, all matches possible
              else:
                  self.matches = self.options[:]
              # No options given, so list files in folder
              if len(self.matches) == 0:
                  self.matches = (glob.glob(text+'*')+[None])
          # return match indexed by state
          try:
              return self.matches[state]
          except IndexError:
              return None

def set_completer(options=[]):
    import readline
    completer = MyCompleter(options)
    readline.set_completer(completer.complete)
    readline.parse_and_bind

