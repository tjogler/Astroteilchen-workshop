#!/usr/bin/env python
import os,sys,glob
import argparse
import utillities as ut
import SEDfitter as sed
import specShapes as spsh
from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np
import testsed
import fitter
import constants as const

'''
def plot(xlist,ylist,error_y=[],error_flaq=[0,0],marker_array=['a','b'],xtitle='E [MeV]',ytitle=r'EF [erg $\mathrm{cm}^{-2} \mathrm{ s}^{-1}$]',lpos='lower left'):
    
    plots the SED with individual colors for each experiment
    plots SED as points with errorbars (set error_flaq =1 for errors)
    plots models as solid curves (when error_flaq=0)
    adds legend with chi2 for models (set marker_array to item name)
    
    line=['-','--',':','-.']
    msymb=['o','s','D','v','^','h','*']
    fig=plt.figure()
    #ax=pyplot.subplots()
    counter=0
    counter_error=0
    counter_line=0
    
    for x in xlist:
        print 'adding plot %s of %s(%s)'%(counter,len(xlist),len(ylist))
        if np.size(x)!=np.size(ylist[counter]):
            print 'Error: x and y not the same dimensions (%i,%i) in point set %i'%(np.size(x),np.size(ylist[counter]),counter)
            exit()
        #print x,ylist[counter]
        if error_flaq[counter]==1:
            mstyle='o'
            #print 'error ' , error_y1
            #plt.loglog(x,ylist[counter],marker=msymb[counter_error],linestyle='None',label=marker_array[counter])
            plt.errorbar(x,ylist[counter],marker=msymb[counter_error],yerr=error_y[counter_error],linestyle='None',label=marker_array[counter])
            #plt.yscale('log')
            #plt.xscale('log')
            counter_error+=1
        else:
            if counter_line >=4:
                counter_line=0
            plt.plot(x,ylist[counter],linewidth=2.5,linestyle=line[counter_line],label=marker_array[counter])
            #plt.yscale('log')
            #plt.xscale('log')
            counter_line+=1
        counter+=1
    #calc min max of the y and x axis)

    counter=0
    counter_error=0
    for er in error_flaq:
        if er ==1:
            ylist.append(ylist[counter]+error_y[counter_error])
            ylist.append(ylist[counter]-error_y[counter_error])
            counter_error+=1
        counter+=1
    
    ymin=ut.min_listofarrays(ylist)*0.8
    ymax=ut.max_listofarrays(ylist)*1.2
    xmin=ut.min_listofarrays(xlist)*0.8
    xmax=ut.max_listofarrays(xlist)*1.2

    #plt.rc('text',usetex=True)
    
    if xmin >0:
        plt.xscale('log')
    if ymin>0:
        plt.yscale('log')
    print 'setting limits x=[%g,%g], y=[%g,%g]'%(xmin,xmax,ymin,ymax)
    plt.ylim(ymin,ymax)
    plt.ylabel(ytitle)
    plt.xlabel(xtitle)
    plt.xlim(xmin,xmax)
    plt.legend(loc=lpos)


'''    

def run_sedfitter(inputfile,collumns,ispec,ipara,iprocess,erange,distance,size,shape,shellth,age,crosspath,testrun,testerror,config,confpath):
    #print crosspath
    print ipara,type(ipara)
    if ipara!=[]:
        ipara=ut.argparse_list_of_list_conv(ipara)
    else:
        fn=open(confpath)
        confDict=eval(fn.read())
        specname=[]
        counter=0
        for p in iprocess:
            specname.append('%s_%s_%i'%(iprocess[counter],ispec[counter],counter))
            counter+=1
        print 'spectra: ',specname
        ipara=ut.conv_dict2list(specname,['a','b','c','d','e','f','g','h'],confDict)
        print 'parameters retrieved from config files:'
        print ipara
        print '\n'
    if erange!=[]:
        erange=ut.argparse_list_of_list_conv(erange)
    else:
        for p in iprocess:
            erange.append([1e6,1e15])


    #print erange,erange[0][0],erange[0][1]

    
    if testrun:
        '''generate test data only works with initial parameter
        provided at the commandline
        '''
        print ispec[0]
        test=testsed.testdata(15,ipara,erange[0][0],erange[0][1],ispec[0],testerror)
        apara=ipara
        print 'par' ,apara
        data=test.make_data()
        #print data
        '''randomize the initial parameters'''
        ipara=np.array(ipara)
        ipara=np.random.normal(ipara,testerror*abs(ipara))
        newerange=[[data[0][0],data[0][np.size(data[0])-1]]]
        print "energyrange changed to data energy range: ",newerange
        Fitroutine=sed.SEDfitter(data,ispec,ipara,newerange,iprocess,sed.source(),crosspath)
    else:
        metadata=[]
        data=[[],[],[],[]]
        ''' read in the data files'''
        if len(inputfile)==1:
            datapart=ut.fill_array_from_txtfile(inputfile[0],collumns," ",8,4.,3,1)
            print type(datapart)
            metadata+=datapart
            data=ut.build_fitdata(data,datapart)
        else:
            for infile in inputfile:
                if inputfile.index(infile)==len(inputfile)-1:
                    datapart=ut.fill_array_from_txtfile(infile,[1,2,3,-1,4,-1]," ")
                    #print 'type: ',type(datapart)
                    #print datapart
                    datapart=ut.conv_magic_sed(datapart)
                    #print len(datapart)
                    #print 'type: ',type(datapart),datapart
                else:
                    datapart=ut.fill_array_from_txtfile(infile,collumns," ",8,2.,3,1)
                    print type(datapart)
                metadata+=datapart
                data=ut.build_fitdata(data,datapart)
   
        #print metadata
        #testspec=sed.spectrum(spsh.test,ipara,data[0][0],data[0][np.size(data[0])-1],np.size(data[0]))
            

    Fitroutine=sed.SEDfitter(data,ispec,ipara,erange,iprocess,sed.source(),crosspath)
    '''SNR specific setting of source parameters'''
    '''W51C'''
    Fitroutine.set_source(distance,size,age,shape,shellth)
    Fitroutine.source.set_radiation_photons(-7.,2.,0.05)
    Fitroutine.set_im(10,0,150)
    Fitroutine.source.set_radiation_densities(0.25,0.84,3.0e-3,0.9)
    
    if config:
        Fitroutine.set_config_dict(confpath)

    modelSpec=Fitroutine.get_model_func()[0]
    modelval=modelSpec(ipara[0])
    
    def pp2T_func(ppFunc,MParticle):
       # newfunc=np.sqrt((ppFunc*const.C)**2+MParticle**2*const.C**4)-MParticle*const.C**2
        newfunc=np.sqrt(ppFunc**2+MParticle**2)#-MParticle
        return newfunc
    
    def dp2dT_func(ppFunc,MParticle,p):
        newfunc=ppFunc/((np.sqrt(p**2+MParticle**2)*p-MParticle**2))
        return newfunc

    
    protonPVal=Fitroutine.parentspec[0].get_xvals()
    protonEnergy=Fitroutine.parentspec[0].get_xvals()
    protonSpecVal=dp2dT_func(Fitroutine.parentspec[0].value(),const.ME_GeV*1.e3,protonPVal)

    finalPar,finalParError,finalChi2,finalVal,func,fineGamma,xFineGamma=Fitroutine.fit(config)
    
    #print 'type func: ',type(func),func
    
    if testrun:
        print 'simulation parameter: ',apara

    '''
    Prepare the plotting of the data and the fit functions
    in a pretty way
    '''

    xvals=[]
    yvals=[]
    yerrors=[]
    yerrflaq=[]
    tags=['Fermi/LAT','MAGIC']
    counter=0
    Fscale=const.MeV2erg
    for p in iprocess:
       ''' if p =="pp_wos":
            tags.append('pp $\chi^2:$ %.1f'%(finalChi2[counter]))
        else: '''   
       tags.append('%s $\chi^2:$ %.1f'%(p,finalChi2[counter]))
       counter+=1
        
    #tags=['Fermi/LAT','MAGIC','%s $\chi^2:$ %.1f'%(finalChi2[0]),'Brems $\chi^2:$ %.1f'%(finalChi2[1]),'IC $\chi^2:$ %.1f'%(finalChi2[2])]
    counter=0
    for i in range(len(metadata)/4):
        xvals.append(metadata[counter])
        yvals.append(metadata[counter+1]*Fscale)
        yerrors.append(metadata[counter+3]*Fscale)
        yerrflaq.append(1)
        counter+=4
    for f in fineGamma:
        if type(f)==list:# required for secondaries in pp process
            for ff in f:
                yvals.append(ff*Fscale)
                tags.append('sub')
                yerrflaq.append(0)
        else:    
            yvals.append(f*Fscale)
            yerrflaq.append(0)
    for xf in xFineGamma:
        if type(xf)==list:# required for secondaries in pp process
            for xxf in xf:
                tags.append('sub')
                xvals.append(xxf)
        else:    
            xvals.append(xf)

    if len(xvals)!= len(yvals):
        print 'ERROR: number of x arrays not equal to y arrays'
        exit()

    tags.append('pp only')
    tags.append('sec brems')
    tags.append('sec IC')
    
    #print xvals,yvals
    #print len(xvals),len(yvals),len(yerrflaq),len(tags)

    fig=ut.plot(xvals,yvals,yerrors,yerrflaq,tags,res=False,onlybands=False,forceLog=True,ytitle=r'E$^{2}$ dN/dE [erg $\mathrm{cm}^{-2} \mathrm{ s}^{-1}$]')
    fig.show()


    ''' 
    xvals=[data[0]]
    finalVal.insert(0,data[1])
            
    for f in fineGamma:
        finalVal.append(f)
        xvals.append(data[0]) # adds the xvals for the rough model fucntions
    for xf in xFineGamma:
        xvals.append(xf)# adds the xvals for the fine model fucntions
    
    if len(finalVal)!=len(xvals):
        print 'ERROR: number of plotted x and y lists is not equal'
        exit()
    
    plot(xvals,finalVal,data[3])
    '''    
    
    #plot([data[0],data[0],func.get_xvals()/1.e6],[data[1],finalVal,func.value()],data[4])
    #plot([data[0],10**np.arange(-2.+3.,4.+3.,0.1),protonEnergy],[data[1],modelval,protonSpecVal],data[4])
    #plot([data[0],data[0],protonEnergy],[data[1],modelval,protonSpecVal],data[4])

    plt.figure()
    #print 'protonPVal: ',protonPVal
    #print 'Ep: ',protonEnergy
    #print 'Fp: ',protonSpecVal
    '''
    print 'test: ',Fitroutine.parentspec[0].get_parameters()
    Fitroutine.parentspec[0].set_precision(np.size(protonEnergy))
    plt.loglog(protonEnergy,Fitroutine.parentspec[0].value())
    plt.loglog(protonEnergy,protonSpecVal)
    plt.show()
    plt.figure()
    plt.loglog(protonEnergy,protonPVal)
    plt.show()
    '''
  
    var=raw_input("Press enter to exit")

parser=argparse.ArgumentParser()
parser.add_argument("--infile",help="file or @filelist containing all input data files",type=str,nargs='+')
parser.add_argument("--columns",help="provide a list for each input file that cspecify the columns to extract the data [x,y,erxup,erxdoen,eryup,erydown] -1 means not there",type=int,nargs='+',default=[1,4,2,-1,5,-1])
parser.add_argument("--spec",help="provide a list for the spectral shape of each process that should be fitted default is pl, possible is pl,bpl,sbpl,plcut",type=str,nargs='+')
parser.add_argument("--par",help="provide a list of initial parameters in the form 'p1 ..pn' for each process to be fitted",type=str,nargs='+',default=[])
parser.add_argument("--proc",help="provide a list of process to be fitted possible: pp,IC,synch,brems",type=str,nargs='+')
parser.add_argument("--erange",help="provide a list of 'emin emax' for each parent particle population (process)",type=str,nargs='+',default=[])
parser.add_argument("--dist",help="distance in pc to the source, default= 2000",type=float,default=2000.)
parser.add_argument("--size",help="provide angluar size of the source in arcsec; default=30",type=float,default=30.)
parser.add_argument("--sourceshape",help="source shape, supported are 'shell' and 'sphere';default = sphere",type=str,default='sphere')
parser.add_argument("--shell",help="provide thickness of shell in pc; default=0.01",type=float,default=0.01)
parser.add_argument("--crpath",help="path where table with pp crossection is saved",type=str,default='./pp_data/')
parser.add_argument("--test",action="store_true",help="this option is for testing purposes only it will disable the data file reading and instead create dummy test data from the spectral shape and parameter you provided, this data will be whiggled with gaussian error of size provided with --testerror; default=False", default=False)
parser.add_argument("--testerror",help="does only work with --test option and is the factor err to calculate the sigma for the gaussian whiggeling sigma=err*dummy_data; default=0.15",type=float,default=0.15)
parser.add_argument("--config",action="store_true",help="if conf option is provided the fit setup is done with configuration files, you can specify the location of your configuration file with --confpath; default=False", default=False)
parser.add_argument("--confpath",help="path+filename to the configuration file, the default is the standard configuration file ./fitConfig.py that can be used as a reference on how to construct the configuration dictionaries",type=str,default='./fitConfig.py')
parser.add_argument("--age",help="age of the source, default 30 ky",type=float,default=3.e4)
args=parser.parse_args()
run_sedfitter(args.infile,args.columns,args.spec,args.par,args.proc,args.erange,args.dist,args.size,args.sourceshape,args.shell,args.age,args.crpath,args.test,args.testerror,args.config,args.confpath)
    
    
