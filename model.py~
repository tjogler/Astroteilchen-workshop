#!/usr/bin/env python
import scipy as sp
import numpy as np
import math
import constants as const
import specShapes as spsh
import crossSection as crs

class model(object):
    '''
    Class that provides the values for each model function
    Contains gamma-spectrum for all availabel processes
    '''

    def __init__(self,spec,process,data,source,ppcross,enconv=1.e3,nucfac=1.845):
        self.spec=spec
        self.proc=process
        self.data=data
        self.source=source
        self.ppcrs=ppcross
        self.conv=1.e3
        self.finalGamma=False
        self.nucFactor=nucfac # enhancement of pp cross due to heavier elements
        
    def set_final_gamma_prec(self,prec):
        '''
        Sets the number of points to evaluate the final gamma model.
        This option is used to get a smooth gamma spectrum evaluated at more than 
        just the data points.
        This option sets also the output of the model function 
        to produce the final gamma spectrum
        '''
        self.gammaPrec=prec
        self.finalGamma=True    
        
    def IC_spec_gamma(self,EdNdE_ph,Eph,EdNdE_e,Ee,en):
       
        dlogE_e = np.log(Ee[1] / Ee[0]) * np.ones_like(Ee)#assumes log spacing of energy bins
        dNe = EdNdE_e * dlogE_e
        #print 'TEST: dlogE_e in IC: ',dlogE_e
        #print 'TEST: Ee[3]/Ee[2]: ',np.log(Ee[3]/Ee[2])
        #print 'TEST: Ee[6]/Ee[5]: ',np.log(Ee[6]/Ee[5])

        dLogE_ph = np.log(Eph[1] / Eph[0]) * np.ones_like(Eph)
        dN_ph = EdNdE_ph * dLogE_ph #changed
        '''print 'dNe: ',dNe
        print 'TEST: Eph: ',Eph
        print 'TEST: Ee: ',Ee
        print 'TEST: en: ',en
        print 'TEST: dLogE_ph: ',dLogE_ph
        print 'TEST: Eph[1]/Eph[0]: ',np.log(Eph[1]/Eph[0])
        print 'TEST: Eph[3]/Eph[2]: ',np.log(Eph[3]/Eph[2])
        print 'TEST: Eph[10]/Eph[9]: ',np.log(Eph[10]/Eph[9])
        print 'TEST: dN_ph: ',dN_ph'''
        
        def EdNdE_gamma(x):
            sigma=crs.sigma_ic(x,Eph,Ee)
            return const.C*const.m2cm*np.dot(dN_ph,np.dot(sigma,dNe))
        
        csec_vec = np.frompyfunc(EdNdE_gamma, 1, 1)
        return csec_vec(en)

    def Brems_spec(self,EdNdE_e,Ee,Eg):
        '''
        Calculates the bremsstrahlungs spectrum for given electron spectrum and ambient density 
        '''
        dlogEe=np.log(Ee[1]/Ee[0])*np.ones_like(Ee)
        dNe=EdNdE_e*dlogEe
        sigma=crs.sigma_brems(Eg,Ee,dNe,self.source.Zeff)
        dNdE_gamma=self.source.nHeff*const.C*const.m2cm*sigma
        return dNdE_gamma

    def Brems_spec_ee(self,EdNdE_e,Ee,Eg):
        ''' 
        Calculates the bremsstrahlungs spectrum for ee bremsstrahlung, assuming a given electron spectrum this formula is only valid E > 5 MeV
        '''
        dlogEe=np.log(Ee[1]/Ee[0])*np.ones_like(Ee)
        dNe=EdNdE_e*dlogEe
        sigma=crs.sigma_brems_ee(Eg,Ee,dNe)
        EdNdE_gamma=self.source.nHeff*const.C*const.m2cm*sigma
        
        return EdNdE_gamma

    def EdQdE_pp(self,dNdp_p, p_p,eg, ID_PARTICLE=0):
        '''
        calculate pp to PARTICLE source function
        INPUT:
        dNdp_p - array_like, shape (n,):
        proton density dN / dp, where p is the momentum (1/GeV/cm^3)
        p_p - array_like, shape (n,):
        proton momenta (GeV)
        n_H - float:
            target gas density (1/cm^3)
        ID_PARTICLE - int:
            particle ID: gamma = 0, electron = 1, positron = 2 etc.
w        OUTPUT:
        E dQ/ dE - function of energy:
            spectrum of produced particles (1/cm^3/s)
        '''
        dNdp_p = np.sqrt(dNdp_p[1:] * dNdp_p[:-1])
        E_p0 = np.sqrt(p_p**2 + const.MP_GeV**2)
        T_p0 = E_p0 - const.MP_GeV
    
    
        E_p = np.sqrt(E_p0[1:] * E_p0[:-1])
        T_p = E_p - const.MP_GeV
        
        dT_p = T_p * np.log(T_p0[1:]/T_p0[:-1])
        
        step = lambda x: (1. + np.sign(x))/2.
        positive = lambda x: x * step(x)
        
        def EdNdE(EE):
            csec = lambda TT: self.ppcrs.pp_dict[ID_PARTICLE](TT, EE/TT) * step(TT - EE)
            csec_vec = np.frompyfunc(csec, 1, 1)
            #print 'step', step(T_p - EE)
            '''print 'p_p: ',p_p
            print 'tp: ',T_p
            print 'dT_p: ',dT_p
            print 'ee: ',EE
            #print 'cs: ',csec(100*EE)
            print 'cs(en): ',csec_vec(T_p)'''
            res = const.C *const.m2cm* self.source.nH * np.sum(csec_vec(T_p) * dNdp_p * dT_p)
            #print 'dNdp_p: ',dNdp_p
            #print 'res: ',res
            return res

        #return EdNdE(eg)
        #print 'dNdp_p: ',dNdp_p
        #print self.ppcrs.pp_dict[0]
        func_vec = np.frompyfunc(EdNdE, 1, 1)
        #print 'output: ', func_vec(eg)
        
        return func_vec(eg)
       # return func_vec 
    

    def pi0_spectrum(self,dNdp_p, p_p,E_g,n_H=1.):
        #print 'dNdp: ',dNdp_p
        #print 'E_g : ',E_g
        kappa_pi = 0.17
        dp_p = p_p[1:] - p_p[:-1]
        dNdp_p = np.sqrt(dNdp_p[1:] * dNdp_p[:-1])

        step = lambda x: (1. + np.sign(x))/2.
        positive = lambda x: x * step(x)

        E_p = np.sqrt(p_p**2 + const.MP_GeV**2)
        E_p = np.sqrt(E_p[1:] * E_p[:-1])
        E_kin = E_p - const.MP_GeV
        kin_mask = step(E_kin**2 / kappa_pi**2 - const.MPi0_GeV**2)

        p_pi = np.sqrt(positive(E_kin**2 / kappa_pi**2 - const.MPi0_GeV**2)) + const.epsilon

        def EdNdE_gamma(E_g):
            E_pi_min = E_g + const.MPi0_GeV**2 / (4 * E_g)
            E_p_min =  const.MP_GeV+ E_pi_min / kappa_pi
            dNpi = const.C*const.m2cm* n_H * crs.sigma_pp(E_kin) * dNdp_p * dp_p
            mask = kin_mask * step(E_p - E_p_min)
            return E_g * 2 * np.sum(dNpi / p_pi * mask)
        EdNdE_gamma_vec = np.frompyfunc(EdNdE_gamma, 1, 1)
        return EdNdE_gamma_vec(E_g)
   
    def pp2T_func(ppFunc,Mparticle):
        newfunc=np.sqrt((ppFunc*const.C)**2+MParticle**2)-MParticle
        return newfunc

    def Tp2pp(self,Tp):
        #converts proton kinectic energy to momentum
        return np.sqrt(Tp*(Tp+2*const.MP_GeV))
        

    def model_val(self,para):
        '''preparing the functor by setting the parameters and energy range for one step of the minimizaton'''
        #print type(self.spec),self.spec[0]
        spectrum=self.spec
        spectrum.set_parameters(para)
        emin=self.spec.logemin
        emax=self.spec.logemax
        #print emin, emax
        #en=np.linspace(emin,emax,spectrum.points)
        
        '''
        Sets the precsion of the parent particle spectrum depending on wheather 
        it is fitting the data or computing the final fine binned parent 
        particle spectrum
        '''
        if self.finalGamma:
            ePoints=np.logspace(np.log10(self.data[0][0]),np.log10(self.data[0][len(self.data[0])-1]),self.gammaPrec)/self.conv
            spectrum.set_precision(1000) #better leave this outside the model so it can be easier adjusted
        else:
            ePoints=self.data[0]/self.conv
            spectrum.set_precision(100) #better leave this outside the model so it can be easier adjusted
        en=np.logspace(emin,emax,spectrum.points)*1.e-9 #in GeV
       # print 'en: ',en
       # print 'N_points: ',len(en)
        EdNdE=spectrum.get_xvals()*spectrum.value()# this is in eV
        #EdNdE2=spectrum.get_xvals()*spectrum.value(en)*1.e-9
        #print 'TEST: EdNdE : ',EdNdE
        #print 'TEST: EdNdE2: ',EdNdE2
        if self.proc=="pp":
            flux_pp=self.nucFactor*self.EdQdE_pp(spectrum.value(),en,ePoints)*ePoints*self.conv
            #secEPoints=np.logspace(np.log10(self.data[0][0]),np.log10(self.data[0][len(self.data[0])-1]),spectrum.points)/1.e3
            spec_secondary=self.nucFactor*const.yr2s*self.source.age*(self.EdQdE_pp(spectrum.value(),en,en,1)+self.EdQdE_pp(spectrum.value(),en,en,2))
            #print 'pp- e+- spec: ',spec_secondary
            flux_sec_brems=ePoints*self.conv*self.Brems_spec(spec_secondary,en,ePoints)
            flux_sec_ic=ePoints*self.conv*self.IC_spec_gamma(self.source.phTot,self.source.Eph,spec_secondary,en,ePoints)
            flux=flux_pp+flux_sec_brems+flux_sec_ic
            #*self.source.sourceVolume*(const.PARSEC*const.m2cm)**3
        elif self.proc=="pp_wos":
            flux=self.nucFactor*self.EdQdE_pp(spectrum.value(),en,ePoints)*self.conv*ePoints
        elif self.proc=="pp_a":
            flux=self.nucFactor*self.pi0_spectrum(EdNdE,en,ePoints,self.source.nH) *ePoints*self.conv
        elif self.proc=="test":
            spectrum.set_precision(np.size(self.data[0]))
            flux=spectrum.value()
        elif self.proc=="ic":
            flux=ePoints*self.conv*self.IC_spec_gamma(self.source.phTot/self.source.Eph,self.source.Eph,EdNdE,en,ePoints)
        elif self.proc=="brems":
            flux=ePoints*self.conv*self.Brems_spec(EdNdE,en,ePoints)
        elif self.proc=="synch":
            flux=-999

        if self.finalGamma:
            '''
            Returns a fine binned Gamma spectrum for the process
            In case of the general pp process it also returns each contributing 
            secondary process as a list [total spectrum, pp, bremsstrahlung, IC] 
            '''
            if self.proc=='pp':
                factor=(1./(4.*np.pi*(self.source.dist*const.PARSEC*const.m2cm)**2.))
                return [flux*factor,flux_pp*factor,flux_sec_brems*factor,flux_sec_ic*factor],ePoints*self.conv
            else:
                return flux*(1./(4.*np.pi*(self.source.dist*const.PARSEC*const.m2cm)**2.)),ePoints*self.conv#*(1./(4.*np.pi*(self.source.dist*const.PARSEC*const.m2cm)**2.))
                
        else:
            return flux/(4.*np.pi*(self.source.dist*const.PARSEC*const.m2cm)**2.)
