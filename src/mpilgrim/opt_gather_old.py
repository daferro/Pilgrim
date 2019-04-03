#!/usr/bin/python2.7
'''
*-----------------------------------*
| Module name:  mpilgrim.opt_gather |
| Last Update:  2019/01/21 (Y/M/D)  |
| Main Author:  David Ferro-Costas  |
*-----------------------------------*
'''

#===============================================================#
import time
import os
import sys
#---------------------------------------------------------------#
import names as PN
import rdwr  as RW
#---------------------------------------------------------------#
import common.Exceptions as     Exc
#---------------------------------------------------------------#
from common.dicts import dpt_im
from   common.files      import write_molden
#from   common.files      import read_gauout
from   common.gauout     import read_gauout
from   common.files      import read_fchk
from   common.files      import read_orca
from   common.files      import read_gtsfile
from   common.files      import write_gtsfile
from ClusterConf import ClusterConf
#---------------------------------------------------------------#
from   common.fncs       import time2human
from   common.fncs       import shift2com
from   common.fncs       import get_symbols
from   common.fncs       import numimag
from   common.fncs       import get_molformula
from   common.fncs       import is_string_valid
from   common.fncs       import time2human
from   common.fncs       import fill_string
#---------------------------------------------------------------#
from   common.Molecule   import Molecule
#---------------------------------------------------------------#
from   common.pgs        import get_pgs
#---------------------------------------------------------------#
from   diverse           import ffchecking
#===============================================================#

#===============================================================#
def known_files(files):
    allowed_ext = ["log","out","fchk","gts"]
    files = [ifile for ifile in files if ifile.split(".")[-1].lower() in allowed_ext]
    return files
#---------------------------------------------------------------#
def data_from_file(filename):
    '''
    Read a file given by the user
    '''
    read_methods = []
    read_methods.append(read_gtsfile) #(a) a gts file
    read_methods.append(read_fchk   ) #(b) a fchk file
    read_methods.append(read_gauout ) #(c) a Gaussian output file
    read_methods.append(read_orca   ) #(d) an Orca output file
    for read_method in read_methods:
        try:
          xcc,atonums,ch,mtp,E,gcc,Fcc,masses,clevel = read_method(filename)[0:9]
          if read_method == read_gtsfile: clevel = ""
          # in case no data in masses
          if masses is None or len(masses) == 0 or sum(masses) == 0.0:
             masses = atonums2masses(atonums)
          return (xcc,atonums,ch,mtp,E,gcc,Fcc,masses,clevel)
        except Exc.FileType:
          continue
        except:
          continue
    return None
#---------------------------------------------------------------#
def analyze_filefolder(ff,dtrack):

    global sthwrong
    global wrong1
    global wrong2

    iblank  = "         "

    # Is ff a folder or a file?
    if os.path.isdir(PN.UFOLDER+ff):
        # Get CTC name & list of files
        ctc_name = ff[:-1]
        files = sorted([ff+filename for filename in os.listdir(PN.UFOLDER+ff)])
        print "         |--> %s"%ff
        iblank += "     "
    else:
        # Get CTC name & list of files
        ctc_name = ff.split(".")[0]
        files = [ff]

    # valid name??
    if not is_string_valid(ctc_name,extra="_"):
        print iblank+"|--> %s (invalid name!)"%ff
        sthwrong = True
        wrong1   = True
        return None

    # only files with known extension
    files = known_files(files)

    # Convert to gts (provisional name)
    count = 1
    gtslist = []
    for ifile in files:
        # read file
        data = data_from_file(PN.UFOLDER+ifile)
        # print info
        if ifile in dtrack.keys():
           print iblank+"|--> %s (in %s)"%(ifile.split("/")[-1],PN.IFILE0)
           continue
        elif data is None:
           print iblank+"|--> %s (reading fails!)"%(ifile.split("/")[-1])
           sthwrong = True
           wrong2   = True
           continue
        else:
           print iblank+"|--> %s (new)"%(ifile.split("/")[-1])
        while True:
              prov_gts = PN.DIR1+"%s.gts_%i"%(ctc_name,count)
              count +=1
              if not os.path.exists(prov_gts): break
        # expand data
        xcc,atonums,ch,mtp,E,gcc,Fcc,masses,clevel = data
        # get point group
        pgroup, rotsigma = get_pgs(atonums,masses,xcc)
        # write gts (provisional)
        if not os.path.exists(PN.DIR1): os.mkdir(PN.DIR1)
        write_gtsfile(xcc,atonums,ch,mtp,E,pgroup,rotsigma,gcc,Fcc,prov_gts,level=clevel)
        gtslist.append( (E,ifile,prov_gts) )

    # sort by energy (only those created)
    gtslist.sort()

    # Rename gts files
    info, idx = [], 1
    for E,ifile,prov_gts in gtslist:
        while True:
              ofile = ctc_name+".%003i.gts"%(idx)
              idx += 1
              if ofile in dtrack.values(): continue
              if not os.path.exists(PN.DIR1+ofile): break 

        os.rename(prov_gts,PN.DIR1+ofile)
        info.append( (E,ifile,ofile) )
    return info
#---------------------------------------------------------------#


#===============================================================#
def ls_struc(dctc):
    '''
    extra option to show what's in .ctc file
    to show: ls ctc
    '''
    ib = "       "
    htable = ["struc name","m.form.","num.ifreqs.","ch","mtp","num.confs."]
    # some numbers and list for nice formatting
    ml1 = max([len(ctc)         for ctc in dctc.keys()  ]+[len(htable[0])])
    ml2 = max([len(CTC._mformu) for CTC in dctc.values()]+[len(htable[1])])
    ll  = [ml1,ml2,len(htable[2]),3,3,len(htable[5])]
    # fill strings in htable
    htable = [fill_string(string,ml) for ml,string in zip(ll,htable)]
    # start loop
    htable  = " "+" | ".join(htable)+" "
    divis   = "-"*len(htable)
    stable  = "\n"
    stable += ib+divis+"\n"
    stable += ib+htable+"\n"
    stable += ib+divis+"\n"
    # separate minima from saddle point and from other (xx)
    ctc_min = sorted([ctc for ctc in dctc.keys() if dctc[ctc]._type == 0])
    ctc_ts  = sorted([ctc for ctc in dctc.keys() if dctc[ctc]._type == 1])
    ctc_xx  = sorted([ctc for ctc in dctc.keys() if dctc[ctc]._type not in [0,1]])
    for ctcs in [ctc_min,ctc_ts,ctc_xx]:
        if ctcs == []: continue
        for ctc in ctcs:
            mformu = dctc[ctc]._mformu
            itcs   = dctc[ctc]._itcs
            ch     = dctc[ctc]._ch
            mtp    = dctc[ctc]._mtp
            sptype = dctc[ctc]._type
            # variable to string
            ncnfs = "%i"%int(sum([weight for itc,weight in itcs]))
            ch    = "%i"%ch
            mtp   = "%i"%mtp
            nif   = "%i"%sptype
            # fill strings
            ltable = [ctc,mformu,nif,ch,mtp,ncnfs]
            ltable = [fill_string(string,ml) for ml,string in zip(ll,ltable)]
            # add to table
            ltable  = " "+" | ".join(ltable)+" "
            stable += ib+ltable+"\n"
        stable += ib+divis+"\n"
    print stable
#===============================================================#


#===============================================================#
# Checkings                                                     #
#===============================================================#
def check_charge(lcharge):
    if len(lcharge) == 1: return lcharge[0]
    else                : return "??"
#--------------------------------------------------------------#
def check_multipl(lmultipl):
    if len(lmultipl) == 1: return lmultipl[0]
    else                 : return "??"
#---------------------------------------------------------------#
def check_nimag(limags):
    if len(limags) == 1: return "%i"%(limags[0])
    else               : return "??"
#---------------------------------------------------------------#
def check_molformula(lmolformu):
    if len(lmolformu) == 1: return lmolformu[0]
    else                  : return "??"
#---------------------------------------------------------------#
def check_gtss(ctc,data_gts):
    ctc_data = data_gts[ctc]
    l_ch     = list(set( [ ch     for (itc,ch,mtp,E,nimag,mformu,pgroup) in ctc_data] ))
    l_mtp    = list(set( [ mtp    for (itc,ch,mtp,E,nimag,mformu,pgroup) in ctc_data] ))
    l_imag   = list(set( [ nimag  for (itc,ch,mtp,E,nimag,mformu,pgroup) in ctc_data] ))
    l_mformu = list(set( [ mformu for (itc,ch,mtp,E,nimag,mformu,pgroup) in ctc_data] ))
    # check elements
    ctc_mformu = check_molformula(l_mformu)
    ctc_type   = check_nimag(     l_imag  )
    ctc_ch     = check_charge(    l_ch    )
    ctc_mtp    = check_multipl(   l_mtp   )
    # inconsistences??
    if "??" not in [ctc_mformu,ctc_type,ctc_ch,ctc_mtp]:
      #print "      num(gts)     = %i"%len(ctc_data)
      #print "      mol.formula  = %s"%ctc_mformu
      #print "      charge       = %i"%ctc_ch
      #print "      multiplicity = %i"%ctc_mtp
      #print "      num(ifreqs)  = %s"%ctc_type
      #print
       isitok = True
    else:
      #print "      num(gts)     = %i"%len(ctc_data)
       print "      ERROR: Inconsistences! Omitting this structure..."
       print "         -------------------------------------------"
       print "           itc  | charge | mtp | n.imag. | m.form.  "
       print "         -------------------------------------------"
       for itc in sorted(ctc_data): # itc = (itc,ch,mtp,E,nimag,mform,pgroup)
           line_data = (itc[0],itc[1],itc[2],itc[4],itc[5])
           print "          %5s |   %2i   | %3i |   %3i   | %-s"%line_data
       print "         -------------------------------------------"
       isitok = False
    return isitok
#---------------------------------------------------------------#
def pgroup2weight(pgroup):
    #if pgroup.lower() == "c1": return 2
    #if pgroup.lower() == "cs": return 1
    return 1
#---------------------------------------------------------------#
def gts_data_and_molden(gtsfile):
    # read gts and prepare Molecule
    molecule = Molecule()
    molecule.set_from_gts(gtsfile)
    molecule.setup()
    # type of PES stationary point?
    nimag  = numimag(molecule._ccfreqs)
    mformu = molecule._mform
    # Get ctc and itc
    ctc,itc,ext = gtsfile.split("/")[-1].split(".")
    # tuple
    gts_tuple = (itc,molecule._ch,molecule._mtp,molecule._V0,nimag,mformu,molecule._pgroup)
    # molden file
    if not os.path.exists(PN.DIR5): os.mkdir(PN.DIR5)
    if is_string_valid(ctc,extra="_"):
       molden  = PN.get_gtsmolden(ctc,itc)
       if not os.path.exists(molden):
           print "      creating %s"%molden
           bmolden = True
           write_molden(molden,molecule._xcc,molecule._symbols,\
                        molecule._ccfreqs,molecule._ccFevecs)
       else:
           bmolden = False
    else:
       molden  = None
       bmolden = False
    return ctc, gts_tuple, molden, bmolden
#===============================================================#


#===============================================================#
def main(idata,status):

    global sthwrong
    global wrong1
    global wrong2
    global wrong3
    global wrong4
    sthwrong    = False
    wrong1      = False
    wrong2      = False
    wrong3      = False
    wrong4      = False

    #---------------------------------------#
    # Read tracking and pif.ctc (if exist)  #
    #---------------------------------------#
    dtrack         ,(fname1,status1) = RW.read_track()
    (dctc,dimasses),(fname2,status2) = RW.read_ctc()
    if status1 == 1: print "  - File '%s' exists and is not empty\n"%fname1
    if status2 == 1: print "  - File '%s' exists and is not empty\n"%fname2

    #---------------------#
    # add data to imasses #
    #---------------------#
    if dimasses == {}: dimasses = dpt_im

    #---------------------------------------#
    # Read user data and generate gts files #
    #---------------------------------------#
    # files and folders
    if os.path.exists(PN.UFOLDER):
       files   = [ff     for ff in os.listdir(PN.UFOLDER) if not os.path.isdir(PN.UFOLDER+ff)]
       folders = [ff+"/" for ff in os.listdir(PN.UFOLDER) if     os.path.isdir(PN.UFOLDER+ff)]
       # only files with known extension
       files = known_files(files)
    else:
       files   = []
       folders = []

    if files+folders != []:
       print "  - Reading user's input data from %s"%PN.UFOLDER
       print
       print "    * Number of files inside %s:"%PN.UFOLDER
       print "      num_files   = %i "%(len(files))
       print "      num_folders = %i "%(len(folders))
       print

       print "      %s"%(PN.UFOLDER)
       for ff in sorted(files)+sorted(folders):
           thedata = analyze_filefolder(ff,dtrack)
           if thedata is None: continue
           for E, ifile, ffile in thedata: dtrack[ifile] = ffile
       print

       #--------------------------------#
       # Write tracking and check files #
       #--------------------------------#
       print "    * Writing/Updating file: %s"%PN.IFILE0
       print
       RW.write_track(dtrack)
       print "    * Checking existency of gts files listed in '%s'"%(PN.IFILE0)
       for ifile, ffile in dtrack.items():
           if not os.path.exists(PN.DIR1+ffile):
             print "      '%s' NOT FOUND! (from '%s')"%(ffile,ifile)
             sthwrong = True
             wrong3   = True
             continue
       print
    else:
       print "  - Folder '%s' does not exist or is empty"%PN.UFOLDER
    print


    #--------------------------------#
    # READ gts files                 #
    #--------------------------------#
    try   : gtsfiles = [PN.DIR1+gts for gts in os.listdir(PN.DIR1) if gts.endswith(".gts")]
    except: gtsfiles = []
    # folder exists?
    if ffchecking([PN.DIR1],[]) == -1 or len(gtsfiles) == 0:
       print "  - Folder '%s' does not exist or is empty"%PN.DIR1
       print "    * Is '%s' empty?"%PN.UFOLDER
       print "    * Does an old '%s' file exist?"%PN.IFILE0
       print
       exit()

    # List gts files
    print "  - Reading gts files in %s"%PN.DIR1
    print
    print "    * num(gts) = %i\n"%len(gtsfiles)

    #--------------------------------------#
    # Analyze them and create molden files #
    #--------------------------------------#
    print "    * molden files will be generated at %s"%PN.DIR5
    data_gts = {}
    num_old, num_new = 0, 0
    for gtsfile in gtsfiles:
        ctc, gts_tuple, molden, bmolden = gts_data_and_molden(gtsfile)
        # check name
        if not is_string_valid(ctc,extra="_"):
            print "      * invalid name (%s)!! Ignoring data..."%ctc
            continue
        # save data
        data_gts[ctc] = data_gts.get(ctc,[]) + [gts_tuple]
        # molden created?
        if bmolden: num_new += 1
        else      : num_old += 1
    print
    print "      # of old molden files: %i"%num_old
    print "      # of new molden files: %i"%num_new
    print 
    print

    #---------------------#
    # Info about each CTC #
    #---------------------#
    print "  - Listing Structures:"
    print
    ml = max([len(ctc) for ctc in data_gts.keys()])
    for ctc in sorted(data_gts.keys()):
        # data for this CTC
        ctc_itcs   = sorted([ (itc,pgroup2weight(pgroup)) \
                        for (itc,ch,mtp,E,nimag,mformu,pgroup)\
                        in   data_gts[ctc]])
        ctc_ch     = data_gts[ctc][0][1]
        ctc_mtp    = data_gts[ctc][0][2]
        ctc_type   = data_gts[ctc][0][4]
        ctc_mformu = data_gts[ctc][0][5]
        ctc_fscal  = 1.0
        ctc_es     = [(ctc_mtp,0.0)]# electronic states
        ctc_dics   = {}
        ctc_diso   = {}
        ctc_anh    = None
        # already exists?
        if ctc in dctc.keys():
           print ("    * %%-%is (already in %%s)"%ml)%(ctc,PN.IFILE1)
           ctc_mformu = dctc[ctc]._mformu
           itcs       = sorted(dctc[ctc]._itcs)
           ctc_ch     = dctc[ctc]._ch
           ctc_mtp    = dctc[ctc]._mtp
           ctc_es     = dctc[ctc]._es
           ctc_type   = dctc[ctc]._type
           ctc_fscal  = dctc[ctc]._fscal
           ctc_dics   = dctc[ctc]._dics
           ctc_anh    = dctc[ctc]._anh
           ctc_diso   = dctc[ctc]._diso
           # to add and to del
           itcs1  = [itc for itc,weight in ctc_itcs]
           itcs2  = [itc for itc,weight in     itcs]
           to_add = [itc for itc in itcs1 if itc not in itcs2]
           to_del = [itc for itc in itcs2 if itc not in itcs1]
           for itc in to_add: print "      - new itc will be added  : %s"%itc
           for itc in to_del: print "      - old itc will be removed: %s"%itc
           # skip if nothing to do
          #if to_add == [] and to_del == []: continue
           # final list
           final_itcs = []
           for itc1,weight1 in ctc_itcs:
               append = True
               for itc2,weight2 in itcs:
                   if itc1 == itc2:
                       final_itcs.append([itc2,weight2])
                       append = False
                       break
               if append: final_itcs.append([itc1,weight1])
           ctc_itcs = final_itcs
        else:
           print ("    * %%-%is"%ml)%ctc
        # checking
        isitok = check_gtss(ctc,data_gts)
        if isitok is False:
           sthwrong = True
           wrong4   = True
           continue
        # save data
        CTC = ClusterConf(ctc)
        CTC._mformu = ctc_mformu
        CTC._itcs   = ctc_itcs
        CTC._ch     = ctc_ch
        CTC._mtp    = ctc_mtp
        CTC._es     = ctc_es
        CTC._type   = ctc_type
        CTC._fscal  = ctc_fscal
        CTC._anh    = ctc_anh
        CTC._diso   = ctc_diso
        CTC._dics   = ctc_dics
        dctc[ctc] = CTC
    print

    # check all ctc has gts files
    for ctc in dctc.keys():
        if ctc in data_gts.keys(): continue
        print "  - No gts file found for '%s'"%ctc
        answer = raw_input("    remove it from %s? (Y/n) "%PN.IFILE1).strip()
        print
        if answer.lower() in ["","y","yes"]: dctc.pop(ctc)

    print "  - Summary table of structures:"
    ls_struc(dctc)

    #----------------#
    # Write ctc file #
    #----------------#
    print "  - Writing/Updating information in %s"%PN.IFILE1
    RW.write_ctc(dctc,dimasses)
    print


    #----------------#
    # EVERYTHING OK? #
    #----------------#
    if sthwrong:
       print "  - ERROR: It seems that something did not go well"
       if wrong1: print "           * invalid name for some file/folder(s)!"
       if wrong2: print "           * reading process failed for a/some file(s) in %s!"%PN.UFOLDER
       if wrong3: print "           * gts file(s) not found!"
       if wrong4: print "           * inconsistences in conformational cluster(s)!"
    else:
       print "  - Everything seems to be OK"
#===============================================================#




