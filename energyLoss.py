#!/usr/bin/env python
import scipy as sp
import numpy as np
import astropy as astro
import math
import constants as const
import specShapes as specsh
import crossSection as crs
import model as mod

from scipy import special

class energyLoss(object):
    ''' class that calculates energy losses for synchrotron, bremsstrhalung
    and IC
    it stores the cooling timescales for each process
    '''
    
    def __init__(self,spec='',proc='',source='',coolingType=''):
        self.spec=spec
        self.proc=proc
        self.source=source
        self.cType=coolingType

    def IC_Edot(self,Eph,Ee,EdNdE_ph):
        f = np.sqrt(Eph[1] / Eph[0])
        dEph = Eph * (f - 1/f)

        E_gb = np.logspace(np.log10(Ee) - 10., np.log10(Ee), 101)
        E_g = np.sqrt(E_gb[1:] * E_gb[:-1])
        dE_g = E_gb[1:] - E_gb[:-1]
        
        sgm = np.zeros_like(Eph)
    
        for i in range(len(E_g)):
            sgm += crs.sigmaIC(E_g[i], Eph, Ee)[:,0] * dE_g[i]


        return const.C*const.m2cm* np.sum(sgm * EdNdE_ph/Eph * dEph)

    def Sync_Edot(self,B,Ee,pitch):
        '''
        calculate the synchrotron energy loss for an electron with energy E_e
        INPUT:
        B - magn field (micro Gauss)
        pitch - sin of alpha (pitch angle)
        E_e - electron energy (GeV)
        OUTPUT:
        Edot - energy loss
        
        '''
        
        def bessel_int_values():
            #dr = 0.00005
            #rs = np.arange(0. + dr/2., 5., dr)
            rsb = np.logspace(-7., 1., 10001)
            rs = np.sqrt(rsb[1:] * rsb[:-1])
            dr = (rsb[1:] - rsb[:-1])
            

            nn = 5./3.
            ks = sp.special.kv(nn, rs)
            k_int = np.zeros_like(ks)

            k_int[0] = ks[0] * dr[0]
            for i in range(1, len(ks)):
                k_int[i] = k_int[i-1] + ks[i] * dr[i]
            k_tot = np.sum(ks * dr)
            k_int = k_tot - k_int

            return rs, k_int

        rs, k_int = bessel_int_values()

        func0 = sp.interpolate.interp1d(rs, k_int, fill_value=0.)
        def func(r):
            if r > 2.e-7 and r < 5.:
                return func0(r)
            else:
                return 0.

        bessel_int = np.frompyfunc(func, 1, 1)

        def synch_norm(B):
            '''
            normalization factor in synchrotron energy loss
            INPUT:
            B - magn field (micro Gauss)
            OUTPUT:
            normalization factor (GeV / s)

            '''
            return np.sqrt(3.) * const.eGauss**3 * 1.e-6 * B / const.ME_GeV * const.erg2GeV**2
        def nu_crit_norm(B):
            norm = 3 * const.H_Erg * const.C*const.m2cm * const.eGauss * 1.e-6 * const.erg2GeV**2* B / (4*np.pi * const.ME_GeV**3 * const.H_eV*1.e-9)
            #norm = 3 * ee * micro * B / (4*np.pi * me**3) * c_light * erg2GeV
            return norm

        def nu_critical(B, sin_al, E):
            '''
            critical frequency
            INPUT:
            B - magnetic field (micro Gauss)
            sin_al = sin(alpha),
            where alpha - angle between magnetic field and electron velocity
            E - electron energy (GeV)
            OUTPUT:
            critical frequency (Hz)
            
            '''
            norm = nu_crit_norm(B)
            
            if type(E) is np.ndarray and type(sin_al) is np.ndarray:
                return norm * np.outer(sin_al, E**2)
            else:
                return  norm * sin_al * E**2
        
        def single_E(E):
            nuc = nu_critical(B, pitch, E)
            #print 'nuc = %.2e GHz' % (nuc / 1.e-9)
            nus = np.logspace(np.log10(nuc)-7., np.log10(nuc)+1., 200)
            f = np.sqrt(nus[1] / nus[0])
            d_nus = nus * (f - 1./f)
            SS = nus / nuc * bessel_int(nus / nuc)
        
            return synch_norm(B) * np.sum(SS * d_nus)
        
        edot_vec= np.frompyfunc(single_E, 1, 1)
        return edot_vec(Ee)
