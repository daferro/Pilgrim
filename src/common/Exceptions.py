#!/usr/bin/python2.7
'''
*----------------------------------*
| Package    :  common             |
| Module     :  Exceptions         |
| Last Update:  2018/10/06 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains classes used to
raise exceptions
'''

class ABORTED(Exception)        : pass
class FileType(Exception)       : pass
class ReadProblem(Exception)    : pass
class UnknownSoft(Exception)    : pass
class ExeNotDef(Exception)      : pass
class ExeNotFound(Exception)    : pass
class CalcFails(Exception)      : pass
class RstDiffTS(Exception)      : pass
class RstDiffPoint(Exception)   : pass
class NoDLEVELdata(Exception)   : pass
class RstDiffVar(Exception)     : pass
class SDpmNoHess(Exception)     : pass
class SDdiffds(Exception)       : pass
class NoICS(Exception)          : pass
class NoTemps(Exception)        : pass
class NoData(Exception)         : pass
class IncompData(Exception)     : pass
class WrongVar(Exception)       : pass
class WrongInIsomass(Exception) : pass
class NoDLEVELpfn(Exception)    : pass
class DLEVELsthWrong(Exception) : pass
class ICFails(Exception)        : pass
class NoTemplateGiven(Exception): pass
class NoReacMol(Exception)      : pass
class LostConformer(Exception)  : pass
class DiffMassConf(Exception)   : pass
class LackOfGts(Exception)      : pass
class NoRateCons(Exception)     : pass
class OnlyMEP(Exception)        : pass
class TSnotFOUND(Exception)     : pass
class NoGTSfile(Exception)      : pass
class FileIsNotGTS(Exception)   : pass




