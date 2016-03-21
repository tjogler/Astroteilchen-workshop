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
     
def sigma_brems(Eg,Ee,Ne,Z=1.):
   


def sigma_brems_ee(Eg,Ee,Ne):
        

def sigma_pp(E_kin):
    
