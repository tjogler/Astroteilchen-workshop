#!/usr/bin/env python
import scipy
from scipy import integrate
import numpy as np
import math
import constants as const
import specShapes as specsh
import crossSection as crss


from scipy import special


class source(object):
    
    ''' class that provides all information abou the source and the 
    ambient enviroment needed for setting up the gamma-ray producing
    processes'''

    _shapes=["sphere","shell"]

    def __init__(self,distance=1000,angularSize=1,age=3.e4,sourceShape='sphere',shellThickness=0):
        '''set source geometry, radiation fields present and  
        ambient medium properties
        The angular Size is given in deg the distance in pc
        The shell thickness is given in pc'''

        self.dist=distance
        self.angSize=angularSize
        self.age=age
        self.sizeParsec=float(2*np.tan(self.angSize*np.pi/180)*self.dist)
        if sourceShape not in self._shapes:
           print "shape %s not yet implemented using sphere instead"%(sourceShape)
           self.sourceShape="sphere"
        else:
           self.sourceShape=sourceShape
        self.rShell=shellThickness
        self._source_volume()
        self.set_am_prop()
        self.ktStar=const.KT_STAR
        self.rhoStar=const.RHO_STAR
        self.ktDust=const.KT_DUST
        self.rhoDust=const.RHO_DUST
        self.ktCMB=const.KT_CMB
        self.rhoCMB=const.RHO_CMB
        self.set_radiation_photons(-6.,2.,0.5)
        self.radiation_fields()

    def _source_volume(self):
        if self.sourceShape=="sphere":
            self.sourceVolume=((self.sizeParsec)**3)/6*np.pi
        elif self.sourceShape=="shell":
            self.sourceVolume=(self.sizeParsec**3-(self.sizeParsec-2*self.rShell)**3)/6*math.pi
        else:
            print "shape not yet implemented using sphere instead"
            self.sourceVolume=(self.angularSize**3)/3*math.pi

    def set_source_geo(self,distance,angularSize,sourceShape='sphere',shellThickness=0):
        self.dist=distance
        self.angSize=angularSize
        if sourceShape not in self._shapes:
           print "shape %s not yet implemented using sphere instead"%(sourceshape)
           self.sourceShape="sphere"
        else:
           self.sourceShape=sourceShape
        self.rShell=shellThickness
        self._source_volume()

    def set_am_prop(self,nH=100.,nIH=80.,Z=1.,T=2.,B=1.):
        '''
        Sets ambient medium properties:
        B is in muGauss
        '''
        self.nH=nH
        self.nIH=nIH
        self.amT=T
        self.Z=Z
        self.B=B
        self.set_eff_prop(self.nH,self.Z)
    
    def set_eff_prop(self,nHeff,Zeff):
        ''' fucntion to set the effective Hydrogen density 
        and effective charge number. 
        Used to account for the abundancy of heavier elements than H'''
        self.nHeff=nHeff
        self.Zeff=Zeff
    
    def set_radiation_densities(self,ktStar,rhoStar,ktDust,rhoDust):
        self.ktStar=ktStar
        self.rhoStar=rhoStar
        self.ktDust=ktDust
        self.rhoDust=rhoDust
        self.radiation_fields()

    
    def set_radiation_photons(self,minEph,maxEph,Eint):
        ''' 
        sets radiation phtoton fields range
        quantities given as log10 values in eV
        evenly spaced in log10 scale
        finally updates all radiationfields with new photon
        energies
        '''
        self.Eph=10**np.arange(minEph,maxEph,Eint)
        self.radiation_fields()


    def radiation_fields(self):
        '''calculates the radiation present in the source
        takes into account CMB,Dust and Starlight'''
      
        #self.Eph=10**np.arange(-7.,2.01,0.05) #radiation field photon energies in eV

        def U_black_body(T):
            '''
            INPUT:
            T - temperature (eV)
            OUTPUT:
            U - energy density of black body radiation (eV/cm^3)
            '''
            aB_eV = 4. * np.pi * 12. * scipy.special.zeta(4. ,1.) / (const.C*const.m2cm * const.H_eV)**3
            return aB_eV * T**4

        def thermal_spectrum(T, U=None):
            '''
            calculate thermal spectrum of photons
            INPUT:
            T - temperature (eV)
            U - energy density of radiation (eV / cm^3) (nu dU / d nu)
            OUTPUT:
            function(Eg) - spectrum of photons as a function of Eg in eV

            '''
            W0 = U_black_body(T) # energy density in the black body with temperature T
            #print 'energy density of CMB ev / cm^3'
            #print W0
            #exit()
            if U is None:
                U = W0
            def func(Eg):
                k = Eg / const.H_eV / (const.C*const.m2cm)
                expf = 1. / (np.exp(Eg/T) - 1.)
                return 8. * np.pi * Eg * k**3 * expf * U / W0
            return func

        self.phStarFunc=thermal_spectrum(self.ktStar,self.rhoStar)
        self.phDustFunc=thermal_spectrum(self.ktDust,self.rhoDust)
        self.phCMBFunc=thermal_spectrum(self.ktCMB,self.rhoCMB)
        self.phStar=self.phStarFunc(self.Eph)
        self.phDust=self.phDustFunc(self.Eph)
        self.phCMB=self.phCMBFunc(self.Eph)
        self.phTot=self.phStar+self.phDust+self.phCMB


        

class spectrum(object):
    
    emin=[]
    emax=[]
    par=[]

    def __init__(self,specfunc,para,emin=1e9,emax=1e15,precision=1000):
        if np.size(specfunc)>1:
            self.val=specfunc[0]
        else:
            self.val=specfunc
        self.par=para
        self.emin=emin
        self.emax=emax
        self.logemin=np.log10(emin)
        self.logemax=np.log10(emax)
        self.points=precision
        self.scaleFact=1.
        if not self._consistancy():
            print 'ERROR: inconsistance in spectrum definition'
            print 'parameters: ',self.par
            

    def _consistancy(self):
        con=True
        if type(self.par)==list:
            for i in self.par:
                if type(i) not in [int,float]:
                    con=False
                    return False
       
        return con


    def get_xvals(self):
        '''returns the x-values as they are used in plotting e10 notation take log10 to obtain the x-values used in the fitting'''
        #x=np.logspace(self.logemin,self.logemax,self.points)
        x=np.logspace(self.logemin,self.logemax,self.points)
        return x

    def get_parameters(self):
        return self.par

    def get_parameter_def(self):
        '''
        returns a list containing a list of the paramerters as first entry
        and a second list explaining if the parameter is linear (0) or 10**p (1)
        '''
        pardef=[]
        pardef.append(self.par)
        pardef.append(specsh.func(self.funcName)[2])
        return pardef

    def integral_whole_range(self,error=0.05):
        '''
        calculates integral using scipy integrate with quad method
        returns value if error on integral less than error
        otherwise it returns -999
        '''
        integr=scipy.integrate.quad(lambda x:self.val(x,self.par),self.emin,self.emax)
        if integr[1]<error*integr[0]:
            return integr[0]
        else:
            print 'ERROR: Integration does not achieve specified uncertainty (%.2f) level'%(error)
            print 'ERROR: I=%.2g eI=%.2g'%(integr[0],integr[1])
            return -999
        
    '''def integral(self,xmin=emin,xmax=emax,ipar=par,step=1.e6):
        print self.par,emin,emax
        xval=np.linspace(xmin,xmax,step)
        dx=abs(xval[0]-xval[1])*np.ones_like(xvals)
        return np.dot(self.value2(xval,ipar),dx)
    '''

    def integral(self,xmin=emin,xmax=emax,ipar=par,error=0.05):
        '''
        calculates integral using scipy integrate with quad method
        returns value if error on integral less than error
        otherwise it returns -999
        '''
        if xmin >100:
            xmin,xmax,ipar[0]=self.scale_erange(xmin,xmax,ipar[0])
        
        integr=scipy.integrate.quad(lambda x:self.val(x,ipar),xmin,xmax)
        
        if integr[1]<error*integr[0]:
            return integr[0]*self.scaleFact
        else:
            print 'ERROR: Integration does not achieve specified uncertainty (%f) level'%(error)
            print 'ERROR: I=%.2g eI=%.2g'%(self.scaleFact*integr[0],self.scaleFact*integr[1])
            return -999
    
    def set_parameters(self,parlist):
        self.par=parlist
    
    def set_erange(self,emin,emax):
        self.emin=emin
        self.emax=emax
        self.logemin=np.log10(emin)
        self.logemax=np.log10(emax)
    
    def set_erange_log(self,emin,emax):
        self.logemin=emin
        self.logemax=emax
        self.emin=10**emin
        self.emax=10**emax

    def set_func_name(self,name):
        self.funcName=name
    
    def set_precision(self,prec):
        self.points=prec
    
    def scale_erange(self,emax,emin,norm):
        self.scaleFact=1./np.log10(emin)
        return emin/emax,emax/emin,norm*self.scaleFact
    
    def value(self,x=np.inf):
        #print 'TEST: type(x) = ',type(x)
        if type(x) != np.ndarray and type(x) != list:
            if x == np.inf:
                '''set spacing of x values as used in the fitting'''
                x=np.logspace(self.logemin,self.logemax,self.points)
                #print 'x in spec.value(): ',x
        return self.val(x,self.par)

    def value2(self,x,par):
        return self.val(x,par)   


