#!/usr/bin/env python
import scipy as sp
import numpy as np
import astropy as astro
import math
import constants as const
import iminuit as minuit 
import specShapes as spsh
import crossSection as crs
import model as mod
import utillities as ut

class fitter(object):
    data=[]
    

    def __init__(self,spec,ipar,process,data,source,confdata,ppcross,method='chi2',verbos=0):
        self.spec=spec
        self.ipar=ipar
        self.proc=process
        self.data=data
        self.source=source
        self.ppcrs=ppcross
        self.conf=confdata
        self.finalGamma=False
        self.model=mod.model
        self.modellist=[]
        self.verbosity=verbos
        self.method=method
        


    def set_final_gamma_prec(self,prec):
        '''Sets the number of points to evaluate the final gamma model. This option is used to get a smooth gamma spectrum evaluated at more than just the data points
        This option sets also the output of the model function to produce the final gamma spectrum'''
        for m in self.modellist:
            m.set_final_gamma_prec(prec)
        #self.gammaPrec=prec
        #self.finalGamma=True
    
    def conv_mval2list(self,mvals):
        '''fucntion that converts the values of the minuit parameters to a list'''
        key=sorted(mvals.keys())
        vpara=[]
        for k in key:
            vpara.append(mvals[k])
        return vpara
    
    def chi2(self,npar):
        if npar <0 or npar >6:
            print 'unsupported number of parameters (%i) for chi^2 function provided'%npar
        if npar==1:
            return self.chi2_1
        elif npar==2:
            return self.chi2_2
        elif npar==3:
            return self.chi2_3
        elif npar==4:
            return self.chi2_4
        elif npar==5:
            return self.chi2_5
        elif npar==6:
            return self.chi2_6
        elif npar==7:
            return self.chi2_7

    def chi2_proto(self,par):
        c2=0.
        c2=(((self.model(par)-self.data[1])**2)/self.data[3]**2)
        c2=c2.sum()
        return c2

    def chi2_1(self,a):
        par=[a]
        c2=0.
        c2=(((self.model(par)-self.data[1])**2)/self.data[3]**2).sum()
        if self.verbosity>0:
            print "chi^2: ",c2
        return c2
    
    def chi2_2(self,a,b):
        par=[a,b]
        c2=0.
        c2=(((self.model(par)-self.data[1])**2)/self.data[3]**2)
        c2=c2.sum()
        if self.verbosiy>0:
            print "chi^2: ",c2
        return c2
    
    def chi2_3(self,a,b,c):
        par=[a,b,c]
        c2=0.
        c2=(((self.model(par)-self.data[1])**2)/self.data[3]**2)
        c2=c2.sum()
        if self.verbosity>0:
            print "chi^2: ",c2
        return c2

    def chi2_4(self,a,b,c,d):
        par=[a,b,c,d]
        c2=0.
        c2=((self.model(par)-self.data[1])**2/(self.data[3]**2))
        #print "type c2: ",type(c2)," c2: ",c2
        c2=c2.sum()
        if self.verbosity>0:
            print "chi^2: ",c2
        return c2
    
    def chi2_5(self,a,b,c,d,e):
        par=[a,b,c,d,e]
        c2=0.
        c2=((self.model(par)-self.data[1])**2/self.data[3]**2).sum()
        if self.verbosity>0:
            print "chi^2: ",c2
        return c2

    def chi2_6(self,a,b,c,d,e,f):
        par=[a,b,c,d,e,f]
        c2=0.
        c2=((self.model(par)-self.data[1])**2/self.data[3]**2).sum()
        if self.verbosity>0:
            print "chi^2: ",c2
        return c2
    
    def chi2_7(self,a,b,c,d,e,f,g):
        par=[a,b,c,d,e,f,g]
        c2=0.
        c2=((self.model(par)-self.data[1])**2/self.data[3]**2).sum()
        if self.verbosity>0:
            print "chi^2: ",c2
        return c2

    def prob_fcn(self,mean,sigma1,sigma2,x):
        '''defines an assymmetric gaussian for
        fitting assymetric errors
        sigma1 for X<mean
        sigma2 for X>mean
        '''
        v1=x<mean
        v2=x>=mean
        val1=1/4./np.pi/sigma1*(np.exp(-((x-mean)/sigma1)**2))*v1
        val2=1/4./np.pi/sigma2*(np.exp(-((x-mean)/sigma2)**2))*v2
        return val1+val2
        
        

    def fit(self,conffile=False,initial=[]):
        newpara=[]
        parerror=[]
        chi2s=[]
        speclist=[]
        modelvallist=[]

        counter=0
        for p in self.proc:
            print "Running process No.%i"%(counter)
            print "starting setting initial parameters for %s spectrum"%(p)
            self.modellist.append(mod.model(self.spec[counter],p,self.data,self.source,self.ppcrs))
            self.model=self.modellist[counter].model_val

            if conffile:
                setup=self.conf['%s_%s_%i'%(p,self.spec[counter].funcName,counter)]
                npar=len(self.ipar[counter])
                counter+=1
                #m=minuit.Minuit(self.chi2(npar),**setup)
                #m.tol=1.e2
            else:
                if initial==[]:
                    initial=self.ipar[counter]
                npar=len(initial)
                print 'number of intial parameters %i'%npar
                setup={}
                counter=0
                parname=['a','b','c','d','e','f','g','h']
                for i in initial:
                    setup[parname[counter]]=i
                    setup['error_%s'%parname[counter]]=0.01
                    limits=ut.make_para_limits(i)
                    setup['limit_%s'%parname[counter]]=limits
                    counter+=1
            m=minuit.Minuit(self.chi2(npar),**setup)
            m.tol=1.e2
            ''' start the fitting process'''    
            print "initial parameters set starting migrad"
            m.migrad()
            # m.hesse()
            newpara.append(self.conv_mval2list(m.values))
            parerror.append(self.conv_mval2list(m.errors))
            chi2s.append(m.fval)
            modelvallist.append(self.model(self.conv_mval2list(m.values)))
            speclist.append(self.spec)
            print 'fitting of process %s done chi^2: %f'%(p,chi2s[len(chi2s)-1])
            

        print 'all %i processes fitted \n '%(counter)
        print '\n\n'
        
        return newpara,parerror,chi2s,modelvallist,speclist


    def set_method(self,method):
        self.method=method
