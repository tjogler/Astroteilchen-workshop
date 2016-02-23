import numpy as np
import constants as const
import numeric as num

class crossSection:
    def __init__(self,filepath,name="pp_tune"):
        self.filepath=filepath
        self.name=name
        self._load_crosssection_table_pptune()

    def _load_crosssection_table_pptune(self):
        csec_dict={}
        Tp=np.load(self.filepath+'Tp.npy')
        IDs = range(7)
        ID_dict = {0:'gamma',
                   1:'elec',
                   2:'posi',
                   3:'nue',
                   4:'numu',
                   5:'antinue',
                   6:'antinumu'
               }

        
        for ID in IDs:
            xx = np.load(self.filepath + 'x_%s.npy' % ID_dict[ID])
            EsigmaE = np.load(self.filepath + 'EsigmaE_%s.npy' % ID_dict[ID])
            #print 'esigmae: ',EsigmaE
            csec_dict[ID] = num.interpolate_linear2d(Tp, xx, EsigmaE)
            #print csec_dict[0](10,.5)
        self.pp_dict=csec_dict
        
        
    def dsigmadx_pg(self,egamma,eproton):
        # this crossection was taken from stefan's script must be validated, returns crossection in m^2
        loge=np.log10(eproton/1.e12) #in log10 TeV
        coeff=(np.ones_like(loge)*0.094)-0.014*loge
        scale=(np.ones_like(loge)*0.062)-0.018*loge
        sigma=(np.ones_like(loge)*42)+8*loge
        dsigmadx=sigma*coeff*(1-egamma)**2.65/(egamma+scale)/egamma
        return dsigmadx*const.MILLIBARN

step = lambda x: (1. + np.sign(x))/2.

def sigma_ic(Eg,Eph,Ee):
    '''
    IC cross section
    Blumenthal and Gould 1970 equation (2.48)
    INPUT:
    Eg - number: final photon energy (GeV)
    Eph - array_like, shape (n,) or a number: radiation field energies (eV)
    Ee - array_like, shape (k,) or a number: electron distribution energies (GeV)
    OUTPUT:
    sigma - array_like, shape (n, k): scattering cross section
    
    '''
    # define useful parameters
    Eph_GeV = Eph / const.GeV2eV
    if not isinstance(Eph_GeV, np.ndarray):
        b = 4 * Eph_GeV * Ee / const.ME_GeV**2
    else:
        b = 4 * np.outer(Eph_GeV, Ee) / const.ME_GeV**2
    z = Eg / Ee
    z = np.minimum(z, 1-const.epsilon)
    x = z / b
    q = x / (1 - z)
    # define the mask to satisfy the conditions for Eg > Eph and Eg < Eg_max = Ee b / (1 + b)
    if isinstance(Eph_GeV, np.ndarray):
        if Eph_GeV[-1] < Eg:
            ph_mask = 1.
        else:
            ph_mask = np.outer(step(Eg - Eph_GeV), np.ones_like(Ee))
    else:
        if Eph_GeV < Eg:
            ph_mask = 1.
        else:
            return 0.

    Eg_max = Ee * b / (1. + b)
    e_mask = step(Eg_max - Eg)
    tot_mask = ph_mask * e_mask
    
    # calculate cross section (Blumenthal and Gould 1970)
    brackets = 2. * q * np.log(q) + (1. + 2. * q) * (1. - q) \
               + 0.5 * z**2 / (1. - z) * (1. - q)
    
    return 3. * const.sigma_Th * x * brackets * tot_mask

def sigma_brems(Eg,Ee,Ne,Z=1.):
    '''
    Bremsstrhalungs cross section according to Blumenthal and Gould 1970 
    Eq 3.1
    Egi in GeV
    Ee in GeV
    '''
    def sigma_ind(Egi):
        maskEef=(Ee-const.ME_GeV)>Egi
        Eef=np.abs(Ee-Egi)+const.epsilon# Eef >0 required because of log using mask to eliminate unphysical contributions
        r02=(const.r0*const.m2cm)**2
        sigma_part1=np.log(2.*Ee*Eef/(Egi*const.ME_GeV))-0.5
        sigma_part2=const.alpha*4.*Z**2*r02/Ee**2*(Ee**2+Eef**2-(2./3.*Ee*Eef))
        sigma=sigma_part1*sigma_part2*maskEef
        #print 'Ef: ',Eef
        #print 'sigma: ',sigma_part1*maskEef
        return np.dot(sigma,Ne)
        

    sigma_vec=np.frompyfunc(sigma_ind,1,1)
    return sigma_vec(Eg)


def sigma_brems_ee(Eg,Ee,Ne):
    ''' bremsstrahlungs cross section on electrons as approximated
    in Baring et al 1999 A1-A4
    it should be only used for Ee>5 MeV below this value a special 
    low energy solution is required
    '''
    def sigma_one(Egi):
        maskEef=(Ee-const.ME)>Egi
        Eef=np.abs(Ee-Egi)+const.epsilon# Eef >0 required because of log using mask to eliminate unphysical contributions
        #print 'Eef: ', Eef
        r02=(const.r0*const.m2cm)**2
        '''
        the first crossection is the same as for eN scattering but Z=1
        '''

        sigma_part1=np.log(2.*Ee*Eef/(Egi*const.ME))-0.5
        #sigma_part1*=step(sigma_part1)*step(Ee-const.ME_GeV-Eg)
        sigma_part2=const.alpha*4.*r02/Ee**2*(Ee**2+Eef**2-(2./3.*Ee*Eef))
        sigma1=sigma_part1*sigma_part2*maskEef
        
        maskEef=(Ee-const.ME)>Egi
        eps_g=Egi/(const.ME)
        gamma_e=(Ee)/const.ME
        '''
        #print 'gamma: ',gamma_e-eps_g,(gamma_e-eps_g)*maskEef
        sigma1_p1=4.*const.alpha*(const.r0*const.m2cm)**2/eps_g
        sigma1_p2=1.+(1./3.-eps_g/gamma_e)*(1.-eps_g/gamma_e)
        sigma1_p3=np.log(2.*gamma_e*np.abs(gamma_e-eps_g+const.epsilon)/eps_g)-0.5#abs insures numerical correct evaluation of log, unphysical solutions are later removed via maskEef
        sigma1=sigma1_p1*sigma1_p2*sigma1_p3*maskEef
        '''
        '''
        Now comes the correction for ee scattering due to recoil of the rest e
        '''
        var1=eps_g<=0.5
        sigma2_p1_1=16.*(1.-eps_g+eps_g**2)*np.log(gamma_e/eps_g)
        sigma2_p1_2=-1./eps_g**2+3./eps_g-4.+4.*eps_g-8.*eps_g**2
        sigma2_p1_3=-2.*(1.-2.*eps_g)*np.log(np.abs(1.-2.*eps_g))
        sigma2_p1_4=(1./4./eps_g**3-1./2./eps_g**2+3./eps_g-2.+4.*eps_g)
        sigma2_p1=sigma2_p1_1+sigma2_p1_2-(sigma2_p1_3*sigma2_p1_4)
        sigma2_p1*=var1

        var2=eps_g>0.5
        sigma2_p2=2./eps_g*((4.-(1./eps_g)+(1./4./eps_g**2))*np.log(2.*gamma_e)-2.+(2./eps_g)-(5./8./eps_g**2))*var2
       # print sigma2_p1
       # print '\n'
       # print sigma2_p2
        sigma2=const.alpha*r02/3./eps_g*(sigma2_p1+sigma2_p2)*maskEef
        #print 'sigma2: ',sigma2
        #print '\n'
        relcorr=1.-(8./3.*(gamma_e-1.)**0.2/(gamma_e+1.)*(eps_g/gamma_e)**(1./3.))
        #print 'corr: ',relcorr,'\n'
        sigmaTot=(sigma2)*relcorr
        
        return np.dot(sigmaTot,Ne)
    sigma_vec=np.frompyfunc(sigma_one,1,1)
    return sigma_vec(Eg)
        

def sigma_pp(E_kin):
        '''
        compute pp cross section (from Aharonian's book)
        INPUT:
        Ekin - kinetic energy of the proton Ekin = E - m (GeV)
        OUTPUT:
        scattering cross section (cm^2)
        '''
        step = lambda x: (1. + np.sign(x))/2.
        positive = lambda x: x * step(x)
        
        E_cut = 1.
        sigma = 30. * (0.95 + 0.06 * np.log(E_kin)) * step(E_kin - E_cut)
        return sigma * const.mb2cm2
