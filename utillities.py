#!/usr/bin/env python

import xml.dom.minidom
import pyfits
import numpy as np
import re
from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch

def get_source_coord_xml(s):
     if s.getAttribute("type").encode()=="PointSource":
         sub=s.getElementsByTagName("spatialModel")[0]
         for p in sub.getElementsByTagName("parameter"):
             if p.getAttribute("name").encode()=="RA":
                 ra=float(p.getAttribute("value").encode())
             elif p.getAttribute("name").encode()=="DEC":
                 dec=float(p.getAttribute("value").encode())
     elif s.getAttribute("type").encode()=="DiffuseSource":
         sub=s.getElementsByTagName("spatialModel")[0]
            #print sub.toxml(),s.toxml()
         spfile=str(sub.getAttribute("file"))
         if not(spfile==""):
             print "name of source map file: ",spfile
             hdu_spfile=pyfits.open(spfile)
             try:
                 ra=float(hdu_spfile[0].header['CRVAL1'])
                 dec=float(hdu_spfile[0].header['CRVAL2'])
                    #print ra,dec
                    #do not add for galactic diffuse model (crude way)
                 if (ra==0 and dec==0):
                     cordsys=hdu_spfile[0].header['CTYPE1']
                     if "GLON" in cordsys:
                         xsize=float(hdu_spfile[0].header['NAXIS1'])*float(hdu_spfile[0].header['CDELT1'])
                         ysize=float(hdu_spfile[0].header['NAXIS2'])*float(hdu_spfile[0].header['CDELT2'])
                         if xsize>=360 or ysize>=180:
                             ra=False
                             dec=False 
             except:
                 ra=False
                 dec=False
         else:
               ra=False
               dec=False

     return ra,dec

def get_sources_par_xml(xmlfile,par):
     xmldoc=xml.dom.minidom.parse(xmlfile)
     sources=xmldoc.getElementsByTagName("source")
     parlist=[]
     NOT=True
     for s in sources:
          sub=s.getElementsByTagName("spectrum")[0]
          for p in par:
               for subsub in sub.getElementsByTagName("parameter"):
                    if subsub.getAttribute("name")==p:
                         parlist.append(subsub.getAttribute("value"))
                         NOT=False
                         break
               if NOT:
                    parlist.append('N/A')
     return parlist

def get_source_par_xml(xmlfile,par,source):
     xmldoc=xml.dom.minidom.parse(xmlfile)
     sources=xmldoc.getElementsByTagName("source")
     parlist=[]
     NOT=True
     for s in sources:
          if source==s.getAttribute("name").encode():
               #print source,s.getAttribute("name").encode()
               sub=s.getElementsByTagName("spectrum")[0]
               for subsub in sub.getElementsByTagName("parameter"):
                    if subsub.getAttribute("name")==par:
                         param=float(subsub.getAttribute("value"))
                         NOT=False
                         break
                    if NOT:
                         param=('N/A')
     return param
               
def import_value(filename,source,quantity):
     para=[]
     if '.dat' in filename:
          print "datfile"
          file=open(filename.strip(),'r')
          res_dict=eval(file.read())
         
        #print quantity
          for q in quantity:
            # might not work with quantity errors and the table related functions
               try:
                    para.append(res_dict[source][q])
               except:
                    para.append('N/A')
                    print 'parameter %s not found for source %s'%(q,source)
        #print para
          file.close()
          if len(para)==1:
               return para[0]
          else:
               return para
     elif '.xml' in filename:
          print 'xmlfile'
          fi=open(filename.strip(),'r+') 
          xmldoc=xml.dom.minidom.parse(fi)
          sources=xmldoc.getElementsByTagName("source")
          for s in sources:
               #print s
               if s.getAttribute("name").encode()==source:
                    for q in quantity:
                         para.append(get_source_par_xml(s,q))
          if len(para)==1:
               return [para[0]]
          else:
               return para
     else:
          print "Error filetype not recognized need .xml or .dat file"
          return 999
          

def get_source_par_xml_list(s,which="all"):
    
    if which not in ['all','free','fixed']:
        print "returning parameter value for %s and if not found N/A"%(which)
    
    NOT=False
    free_para=[]
    fix_para=[]
    sub=s.getElementsByTagName("spectrum")[0]
    for p in sub.getElementsByTagName("parameter"):
         if p.getAttribute("name").encode()==which:
              para_val=p.getAttribute("value").encode()
              NOT=True
         if p.getAttribute("free").encode()=='1':
              free_para.append( p.getAttribute("name"))
         else:
              fix_para.append( p.getAttribute("name"))
    if which=='all':
         return free_para,fix_para
    elif which=='free':
         return free_para
    elif which=='fixed':
         return fix_para
    elif NOT:
         return para_val
    else:
        print 'ERROR in function get_source_par_xml(s,which="all")'
        return -999
    
def getSpecParXML(source):
    para=[]
    spec=source.getElementsByTagName("spectrum")[0] 
    for s in spec.getElementsByTagName("parameter"):
        st=''
        st=s.getAttribute("name").encode()
        st+='='
        st+=s.getAttribute("value").encode()
        st+=' +/- '
        st+=s.getAttribute("error").encode()
        st+=' ,free='
        st+=s.getAttribute("free").encode()
        para.append(st)
    return para


def get_path_from_filestr(filestring):
     fparts=filestring.split('/')
     filename = fparts[len(fparts)-1]
     path=filestring[:len(filestring)-len(filename)]

     return path


def calc_syst(ref,others):
     fsum=0
     for o in others:
          fsum+=(ref-o)**2
     fsum=fsum/len(others)
     return np.sqrt(fsum)


def fill_array_from_txtfile(fname,col=[1,2,-1,-1,-1,-1],sep=" ",cutcol=-1,cutval=-999,header=0,trailer=0):
     '''
     reads in files and returns a list of numpy arrays of 
     the form [x,y,xerup,xerdown,yerup,yerdown] 
     were x and y are followed by any positive error column flag see below
      -1 or any negative number ignores this collumn,
     used if only y errors are provided or no errors at all
     '''
     cols=[]
     for cc in col:
          cols.append(int(cc))
          
     if len(cols) <6:
          for i in range(len(cols),6):
               cols.append(-1)
          
          
     infile=open(fname,'r+')
     if trailer <1:
          lines = infile.readlines()[header:]
     else:
          lines = infile.readlines()[header:len(infile.readlines())-trailer]
     if len(lines)<1:
          print "no lines to read"
          return -999
     #change indices for x and y to numpy counting starting with 0
     nitems=len(lines)
     results=[]
     res=[]
     counter=0
     cols2=cols[:]
     for c in cols2:
          if c >-1:
               results.append(np.zeros(nitems))
               cols2[counter]-=1
          counter+=1

     """ if cols[2] >0:
          xercol-=1
          xerval=np.zeros(nitems)
     if yercol >0:
          yercol-=1
          yerval=np.zeros(nitems)
     """   
     counter=0
     for l in lines:
          part=re.split(sep,l)
          part = [x for x in part if x != ''] #added to remove blancs when using regular expressions should be properly fixed one time
          used=0
          if (cutcol>-1) and (float(part[cutcol-1])>cutval):
               for c in cols2:
                    if c >-1:
                         results[used][counter]=float(part[c])
                         used+=1
               counter+=1
          elif cutcol==-1:
               for c in cols2:
                    if c >-1:
                         #print part[c]
                         results[used][counter]=float(part[c])
                         used+=1
               counter+=1
    
     #Now remove all unused array positions due to the cut criteria
     for r in results:
          res.append(r[0:counter])
     infile.close()
    
     return res


def max_listofarrays(flist):
     maximum=[]
     for f in flist:
          maximum.append(max(f))

     return max(maximum)

def min_listofarrays(flist):
     minimum=[]
     for f in flist:
          minimum.append(min(f))

     return min(minimum)

def argparse_list_of_list_conv(inlist,delimiter=' '):
     '''
     function that conversts each string of a list of strings
     into a list of floats
     this is required to parse multiple arguments as multiple lists 
     from the commandline with argparse
     the function returns a list conataining the elements of each string
     as a sub list
     '''

     result=[]
     for item in inlist:
          sublist=[]
          for subitem in item.split(delimiter):
               sublist.append(float(subitem))
          result.append(sublist)

     return result

def conv_magic_sed(mdata):
     mdata[0]=mdata[0]*1.e3  # converts energy to MeV
     mdata[2]=mdata[2]*1.e3  # converts energy to MeV
     mdata[1]=mdata[1]*(mdata[0]/1.e6)*mdata[0] # converts flux to E^2F in MeV
     mdata[3]=mdata[3]*(mdata[0]/1.e6)*mdata[0]
     return mdata

def conv_hess_sed(mdata,cons=True):
     '''
     converts the assymetric errors to symmetric ones 
     using the smaller error if conserv is set to True otherwise
     uses the arithmetric mean
     output in MeV
     '''
     e1=abs(mdata[1]-mdata[2])
     e2=abs(mdata[1]-mdata[3])
     if cons:
          er1= e1<e2
          er2= e1>=e2
          e=e1*er1+er2*e2
     else:
          e=(e1+e2)/2

     mdata[1]=mdata[1]*mdata[0]**2*1.e6
     mdata[2]=e*mdata[0]**2*1.e6
     mdata[3]=e*mdata[0]**2*1.e6
     mdata[0]=mdata[0]*1.e6
     return mdata

def conv_marianne_data(mdata,cons=True):
     '''
     converts the assymetric errors to symmetric ones 
     using the smaller error if conserv is set to True otherwise
     uses the arithmetric mean
     output in MeV
     '''
    
     e1=abs(mdata[1]-mdata[2])
     e2=abs(mdata[1]-mdata[3])
     if cons:
          er1= e1<e2
          er2= e1>=e2
          e=e1*er1+er2*e2
     else:
          e=(e1+e2)/2
         
     mdata[2]=e
     mdata[0]*=1.e3
     return mdata     


def build_fitdata(data,datapart):
     result=[]
     counter=0
     for d in data:
          res=np.hstack((d,datapart[counter]))
          res.ravel()
          counter+=1
          result.append(res)
     return result


def make_para_limits(para):
     ''' 
     function that calculates useful limit for a fit routine
     parameter limits will be +/- and about 5X large than the initial parameter 
     if the parameter is < 0 the value will be +/- 4
     '''
     limits=np.arange(2)
     if abs(para) <1:
          limits[0]=-4
          limits[1]=4
     elif para >=1:
          limits[0]=-5.*para
          limits[1]=5.*para
     elif para <= -1.:
          limits[0]=5.*para
          limits[1]=-5.*para
          
     return limits 

def count_elements(vlist):
     '''
     function that counts how often each element is in a list and save the
     information in another dimension
     implemented for 1 dimension right now
     '''
     countlist=[]
     
     for l in vlist:
          countlist.append(vlist.count(l))

     return [vlist,countlist]
          

def get_para_from_dict(dic,keynam):
          key=dic.keys()
          vpara=list(np.zeros(len(keynam)))
          countElements=0
          for k in key:
               counter=0
               for kn in keynam:
                    if k==kn:
                         vpara[counter]=dic[k]
                         countElements+=1 #needed to know how many keys are in keyname
                    counter+=1
          return vpara[:countElements] # returns a list with only contained keys that can be less than the number provided in keynam

def conv_dict2list(dictnames,keynames,metadict):
     '''fucntion that converts the values of the dictonary entries
     provided in keynames into a list
     in case that a metadict is provided it makes a list for each dictname
     then it places each of these lists as an element in a metalist and 
     returns the  metalist
     '''
     
     totallist=list(np.ones(len(dictnames)))
     if type(metadict)==dict:
          dictaskeys=metadict.keys()
          #print dictaskeys
          for di in dictaskeys:
               counter=0
               #print 'di: ',di
               #print dictnames
               for name in dictnames:
                    if di==name:
                         #print 'name',name
                         #print 'counter: ',counter
                         totallist[counter]=get_para_from_dict(metadict[di],keynames)
                    counter+=1
               
          return totallist

     else:
          return get_para_from_dict(dictnames,keynames)

def fill_between(x, y1, y2=0, ax=None, **kwargs):
    """Plot filled region between `y1` and `y2`.

    This function works exactly the same as matplotlib's fill_between, except
    that it also plots a proxy artist (specifically, a rectangle of 0 size)
    so that it can be added it appears on a legend.
    """
    ax = ax if ax is not None else plt.gca()
    ax.fill_between(x, y1, y2, **kwargs)
    p = plt.Rectangle((0, 0), 0, 0, **kwargs)
    ax.add_patch(p)
    return p

def replace_val(listoflist=[],repVal=1.e-32,val=0.):
     if type(listoflist) != list:
          print "ERROR: No list of list provided: ",listoflist
     for l in listoflist:
          if type(l) != np.ndarray:
               print "ERROR: No array in list: ",l
          l[l<=val]=repVal
     return listoflist

def plot(xlist,ylist,error_y=[],error_flaq=[0,0],marker_array=['a','b'],bands=[],combine=[],xtitle='E [MeV]',ytitle='EF',lpos='lower left',res=True,onlybands=True,forceLog=False,yscale=[1.e-16,1.e-9]):
     '''
     plots the SED with individual colors for each experiment
     plots SED as points with errorbars (set error_flaq =1 for errors)
     plots models as solid curves (when error_flaq=0)
     
     adds legend with chi2 for models (set marker_array to item name)
     '''
     line=['-','--',':','-.']
     msymb=['o','s','D','v','^','h','*']
     fig=plt.figure()
     
    
     if res:
          left,width=0.2,0.8
          rect1 = [left, 0.3, width, 0.6]
          rect2 = [left, 0.3, width, 0.2]
          ax1 = fig.add_axes(rect1)  #left, bottom, width, height
          ax2 = fig.add_axes(rect2, sharex=ax1)
     else:
          ax1 = fig.add_axes([0.175,0.15,0.8,0.8])  #left, bottom, width, height
     
     counter=0
     counter_error=0
     counter_line=0
     
     '''
     prepare sublist that should be shown as band
     '''

     if combine!=[]:
          checkXrange=xlist[combine[0]:combine[1]+1]
          cx=0
          for x in checkXrange:
               if all(x==checkXrange[cx]):
                    cx+=1
               else:
                    print "ERROR: combining different xranges not possible in band plotting option"
                    print "ERROR: x1=%s"%x
                    print "ERROR: x2=%s"%checkXrange[cx]
                    exit()
     
          ce=0
          ycomb=ylist[combine[0]:combine[1]+1]
          ycombmin=np.ndarray(1)# cannot create empty array must remove first element later
          #print ycombmin
          ycombmax=np.ndarray(1)
          if any(error_flaq[combine[0]:combine[1]+1]):
               print 'errors found'
               yerror=error_y[combine[0]:combine[1]+1]
               elist=error_flaq[combine[0]:combine[1]+1]
               for e in elist:
                    #print yerror
                    if e==1:
                         ycombmin=np.append(ycombmin,ycomb[ce]-yerror[ce])
                         ycombmax=np.append(ycombmax,ycomb[ce]+yerror[ce])
                         ce+=1
                    else:
                         ycombmin=np.append(ycombmin,ycomb[ce])
                         ycombmax=np.append(ycombmax,ycomb[ce])
                         ce+=1
          else:
               print 'no errors found'
               ycombmin=np.append(ycombmin,ycomb[ce])
               ycombmax=np.append(ycombmax,ycomb[ce])
               ce+=1
          #print 'min: ',ycombmin
          #print 'max: ',ycombmax
          ycombmin=ycombmin[1:]
          ycombmax=ycombmax[1:]
          #print ycombmin,len(ycomb[0]),ce

          minlist=ycombmin.reshape(ce,len(ycomb[0])).transpose()
          maxlist=ycombmax.reshape(ce,len(ycomb[0])).transpose()

          

          lowerbound=np.ndarray(1)
          upperbound=np.ndarray(1)
          #print 'min: ',minlist
          for a in minlist:
               #print 'lower: ',lowerbound
               lowerbound=np.append(lowerbound,min(a))
          #print 'max: ',maxlist
          for b in maxlist:
               #print 'upper: ',upperbound
               upperbound=np.append(upperbound,max(b))
          lowerbound=lowerbound[1:] # removing the first element that does not belong to list
          upperbound=upperbound[1:]
          bands.append([lowerbound,upperbound,xlist[combine[0]]])
          #print bands              
          if onlybands:
               if combine[0]==0:
                    combine[0]=1
               xlist=xlist[:combine[0]-1]+xlist[combine[1]+1:]
               ylist=ylist[:combine[0]-1]+ylist[combine[1]+1:]
               error_y=error_y[:combine[0]-1]+error_y[combine[1]+1:]
               error_flaq=error_flaq[:combine[0]-1]+error_flaq[combine[1]+1:]
     #print 'xlist: ', xlist
     #print 'ylist: ', ylist
     #print 'elist: ',error_y

    # for x,y,e in xlist,ylist,error_y:
     #     print 'lenght: ',len(x),len(y),len(e)

     if forceLog:
          replace_val(ylist,yscale[0])
          replace_val(error_y,yscale[0])

     for x in xlist:
         
          print 'adding plot %s of %s(%s)'%(counter,len(xlist),len(ylist))
          if np.size(x)!=np.size(ylist[counter]):
               print 'Error: x and y not the same dimensions (%i,%i) in point set %i'%(np.size(x),np.size(ylist[counter]),counter)
               exit()
        
          if error_flaq[counter]==1:
               mstyle='o'
               if len(ylist)<=counter-1:
                    print 'ylist to short %s %i'%(ylist,counter)
               if len(msymb)<=counter_error-1:
                    print 'msymb to short %s %i'%(msymb,counter_error)
               if len(marker_array)<=counter-1:
                    print 'marker_array to short %s %i'%(marker_array,counter)
               if len(error_y)<=counter_error-1:
                    print 'error_y to short %s %i'%(error_y,counter_error)
               if len(x)<=counter-1:
                    print 'x to short %s %i'%(x,counter) 
               ax1.errorbar(x,ylist[counter],marker=msymb[counter_error],yerr=error_y[counter_error],linestyle='None',label=marker_array[counter])
               counter_error+=1
          else:
               if counter_line >=4:
                    counter_line=0
               #print x, ylist[counter],counter,counter_line,line[counter_line],marker_array[counter]
               ax1.plot(x,ylist[counter],linewidth=2.5,linestyle=line[counter_line],label=marker_array[counter])
               if res:
                    ax2.plot(x,abs(ylist[counter]-ylist[0])/ylist[0],linewidth=2.5,linestyle=line[counter_line],label=marker_array[counter])
               counter_line+=1
          counter+=1
     
     counter_error=0
     count=0
     for er in error_flaq:
          if er ==1 and len(ylist)>0:
               ylist.append(ylist[count]+error_y[counter_error])
               ylist.append(ylist[count]-error_y[counter_error])
               counter_error+=1
          count+=1  
     
     #print 'ylist: ',ylist
     

     for b in bands:
          print 'plotting %i bands to plot'%len(bands)
          #print len(b[2]),len(b[0]),len(b[1])
          #print b[2],b[0],b[1]
          if forceLog:
               replace_val(b,yscale[0])
          fill_between(b[2], b[0], b[1],ax=ax1,alpha=0.3,antialiased=True,label=marker_array[counter])
          xlist.append(b[2])
          ylist.append(b[0])#required to get the scale setting right
          ylist.append(b[1])#required to get the scale setting right
          counter+=1
          #facecolor='#089FFF'          

     

     ymin=min_listofarrays(ylist)*0.8
     ymax=max_listofarrays(ylist)*1.2
     xmin=min_listofarrays(xlist)*0.8
     xmax=max_listofarrays(xlist)*1.2

     if xmin >0:
          ax1.set_xscale('log')
          if res:
               ax2.set_xscale('log')
     if ymin>0:
          ax1.set_yscale('log')
     if forceLog:
          replace_val(ylist,1.e500,yscale[0])
          ymin=min_listofarrays(ylist)*0.8
          print 'some errors go to negative values, truncate y-axis'
          ax1.set_yscale('log')
          
     
     print 'setting limits x=[%f,%f], y=[%e,%e]'%(xmin,xmax,ymin,ymax)
     ax1.set_ylim(ymin,ymax)
     ax1.set_ylabel(ytitle,fontsize=20)
     if res:
          ax2.set_xlabel(xtitle)
          ax2.set_ylabel(r'$\Delta\,$ '+ytitle)
     ax1.set_xlim(xmin,xmax)
     ax1.legend(loc=lpos)
     ax1.set_xlabel(xtitle,fontsize=20)
     ax1.tick_params(axis='x', labelsize=18)
     ax1.tick_params(axis='y', labelsize=18)
     ax1.tick_params(which='major',size=10)
     ax1.tick_params(which='minor',size=5)
     return fig


def convert_energy_units(inunits,outunits,value):
     '''
     function to convert energy units
     it returns the converted value and the conversion factor
     '''
     if inunits== outunits:
          return value,1

     units=['eV','keV','MeV','GeV','TeV','erg']

     if inunits in units:
          i=units.index(inunits)
          #print 'i: ',i
     else:
          print "ERROR: input units not defined"
          exit()
     if outunits in units:
          j=units.index(outunits)
          #print 'j: ',j
     else: 
          print "ERROR: output units not defined"
          exit()

     table=[[1.,1.e-3,1.e-6,1.e-9,1.e-12,6.24e-11],[1.e3,1.,1.e-3,1.e-6,1.e-9,6.24e-8],[1.e6,1.e3,1.,1.e-3,1.e-6,6.24e-5],[1.e9,1.e6,1.e3,1.,1.e-3,6.24e-2],[1.e12,1.e9,1.e6,1.e3,1.,6.24e1],[6.24e11,6.24e8,6.24e5,6.24e2,6.24e1,1.]]
     
     return value*table[i][j],table[i][j]
     
