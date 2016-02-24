#!/usr/bin/env python
import os,sys,glob
import argparse
import utillities as ut
import specShapes as spsh
import sedObj as sed
from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np
import constants as const
import model


def plot_spectrum(ispec,ipara,iprocess,erange,distance,size,shape,shellth,age,crosspath):
    #print ipara,type(ipara)
    if ipara!=[]:
        ipara=ut.argparse_list_of_list_conv(ipara)[0]
    specfunc=spsh.func(ispec[0])[0]
    #print specfunc(np.linspace(erange[0],erange[1],10),ipara)
    '''SNR specific setting of source parameters'''
    '''W51C'''
    Source=sed.source(distance,size,age,shape,shellth)
    Source.set_radiation_photons(-7.,2.,0.05)#energy range and binning for target photons
    Source.set_am_prop(nH=10.,nIH=0.,Z=1.) # inter stellar medium parameters
    Source.set_radiation_densities(0.25,0.84,3.0e-3,0.9) # radiation field desities and temperatures for IR and UV photonfields (CMB is fixed and cannot be changed)

    spec=sed.spectrum(specfunc,ipara,erange[0],erange[1],precision=100)
    #print spec
    gammaEmission=model.model(spec,iprocess[0],Source,crosspath)
    gammaEmission.set_final_gamma_prec(1000)
    fineGamma,xFineGamma=gammaEmission.model_val(ipara)
    
    #print 'type func: ',type(func),func
    
    
    '''
    Prepare the plotting of the functions
    in a pretty way
    '''

    xvals=[]
    yvals=[]
    yerrors=[]
    yerrflaq=[]
    tags=['bla','blub']
    counter=0
    Fscale=const.MeV2erg
    for p in iprocess:
       tags.append('%s'%(p))
       counter+=1
        
    counter=0
    for f in fineGamma:
        if type(f)==list:# required for secondaries in pp process
            for ff in f:
                yvals.append(ff*Fscale)
                yerrflaq.append(0)
        else:    
            yvals.append(f*Fscale)
            yerrflaq.append(0)
        
    for xf in xFineGamma:
        if type(xf)==list:# required for secondaries in pp process
            for xxf in xf:
                xvals.append(xxf)
        else:    
            xvals.append(xf)

    if len(xvals)!= len(yvals):
        print 'ERROR: number of x arrays not equal to y arrays'
        exit()

    
    #print 'x: %s'%xvals
    #print 'y: %s'%yvals
        
    fig=ut.plot([xvals],[yvals],error_flaq=yerrflaq,marker_array=tags,res=False,onlybands=False,forceLog=False,ytitle=r'E$^{2}$ dN/dE [erg $\mathrm{cm}^{-2} \mathrm{ s}^{-1}$]')
    fig.show()


    #plt.figure()
  
    var=raw_input("Press enter to exit")

parser=argparse.ArgumentParser()
parser.add_argument("--spec",help="provide a list for the spectral shape of each process that should be plotted default is pl, possible is pl,bpl,sbpl,plcut",type=str,nargs='+')
parser.add_argument("--par",help="provide a list of initial parameters in the form 'p1 ..pn' for each process to be plotted",type=str,nargs='+',default=[])
parser.add_argument("--proc",help="provide a list of process to be plotted possible: pp,IC,synch,brems",type=str,nargs='+')
parser.add_argument("--erange",help="provide a list of 'emin emax' for each parent particle population (process)",type=float,nargs='+',default=[])
parser.add_argument("--dist",help="distance in pc to the source, default= 2000",type=float,default=2000.)
parser.add_argument("--size",help="provide angluar size of the source in arcsec; default=30",type=float,default=30.)
parser.add_argument("--sourceshape",help="source shape, supported are 'shell' and 'sphere';default = sphere",type=str,default='sphere')
parser.add_argument("--shell",help="provide thickness of shell in pc; default=0.01",type=float,default=0.01)
parser.add_argument("--crpath",help="path where table with pp crossection is saved",type=str,default='./pp_data/')
parser.add_argument("--age",help="age of the source, default 30 ky",type=float,default=3.e4)
args=parser.parse_args()
plot_spectrum(args.spec,args.par,args.proc,args.erange,args.dist,args.size,args.sourceshape,args.shell,args.age,args.crpath)
    
    
