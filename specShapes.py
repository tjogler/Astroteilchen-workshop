import numpy as np

def func(name):
    if name=="pl":
        return [pl,3,[1,0,1]]
    elif name=="bpl":
        return [bpl,4,[1,0,0,1]]
    elif name=="plcut":
        return [plcut,4,[1,0,1,1]]
    elif name=="bplcut":
        return [bplcut,5,[1,0,0,1,1]]
    elif name=="sbpl":
        return [sbpl,5,[1,1,0,1,0]]
    elif name=="sbplcut_fail":
        return [sbplcut_fail,6,[1,1,0,1,0,1]]
    elif name=="sbplcuty":
        return [sbplcuty,6,[1,1,0,1,0,1]]
    elif name=="sbplcut":
        return [sbplcut,7,[1,1,0,1,0,1,1]]
    elif name=="test":
        return [test,2,[1,1]]

def _integral(xvals,func,para,binning=1000.):
    '''
    calculates the integral of a function.
    The xrange is recomputed using 1000 bins
    it is assumed that the normalization of func is
    provided in the para[0]
    x is assumed to be log10 spaced
    '''
    if len(xvals)<2:
        print 'integral requires bounds but no x range supplied'
        exit()
    xF=np.linspace(xvals[0],xvals[len(xvals)-1],binning)
    #dx=abs(np.log10(xF[0]/xF[1]))*np.ones_like(xF)
    dx=abs(xF[0]-xF[1])#*np.ones_like(xF)
    para[0]=0.  #set any normalization parameter to 1
    yvals=func(xF,para)
    #integral=np.dot(yvals,dx)
    integral=yvals.sum()*dx
    return integral

def pl(x,parlist):
    if len(parlist)!=3:
        print "function powerlaw needs 3 parameters but %f given"%len(parlist)
    return (10**parlist[0])*(x/10**parlist[2])**parlist[1]

def bpl(x,parlist):
    if len(parlist)!=4:
        print "function broken powerlaw needs 4 parameters but %f given"%len(parlist)
    v1= x< float(10**parlist[3])
    v2= x>=float(10**parlist[3])
    val1= 10**(parlist[0])*(x/10**parlist[3])**parlist[1]*v1
    val2= 10**(parlist[0])*(x/10**parlist[3])**parlist[2]*v2
    return val1+val2

def plcut(x,par):
    if len(par)<4:
        print "function powerlaw with exponential cut-off needs at least  4 parameters but %f given"%len(par)
    #print 'x as in plcut: ',x
    val2=10**par[0]*(x/10**par[1])**par[2]*np.exp(-x/max(10**par[3],1))
    #val2=(par[0]-(par[2]*(x-par[1]))-10**(x/par[3])*np.log10(np.exp(1)))
    #print 'val2: ',val2
    return val2

def bplcut(x,par):
    if len(par)<5:
        print 'function broken power law with exponential cut-off needs 5 parameters but %i given'%len(par)
    v1= x< float(10**par[3])
    v2= x>=float(10**par[3])
    val1= 10**(par[0])*(x/10**par[3])**par[1]*v1
    val2= 10**(par[0])*(x/10**par[3])**par[2]*np.exp(-x/max(10**par[4],1))*v2
    return val1+val2

def sbpl(x,par):
    if len(par)<5:
        print ' function smooth-broken power law needs 5 parameters but %i given'%len(par)
    
    val1=10**par[0]*(x/10**par[1])**par[2]
    val2=(1.+(x/10**par[3])**par[4])**(-1.)
    
    return val1*val2


def sbplcut_fail(x,par):
    if len(par)<6:
        print ' function smooth-broken power law with exponential cut-off needs 6 parameters but %i given'%len(par)
    
    val1=10**par[0]*(x/10**par[1])**par[2]
    val2=(1.+(x/10**par[3])**par[4])**(-1.)
    val3=np.exp(-x/max(10**par[5],1)) #exponential cut-off
    return val1*val2*val3

def sbplcut(x,par):
    if len(par)<7:
        print ' function smooth-broken power law with exponential cut-off needs 7 parameters but %i given'%len(par)
    
    val1=10**par[0]*(x/10**par[1])**par[2]
    val2=(1.+(x/10**par[3])**((par[2]-par[4])/par[6]))**(-par[6])
    val3=np.exp(-x/max(10**par[5],1)) #exponential cut-off
    return val1*val2*val3

def sbplcuty(x,par):
    if len(par)<6:
        print ' function smoot-broken power law with exponential cut-off needs 6 parameters but %i given'%len(par)
    
    val1=10**par[0]*(x/10**par[1])**par[2]
    val2=(1.+(x/10**par[3])**2)**(-par[4]/2.)
    val3=np.exp(-x/max(10**par[5],1)) #exponential cut-off
    return val1*val2*val3

'''def plcut(x,parlist=[1,1,1,1,1,1]):
    if len(parlist)<4:
        print "function powerlaw with exponential cut-off needs at least  4 parameters but %f given"%len(parlist)
    elif len(parlist)>6:
        print "function powerlaw with exponential cut-off takes not more than 6 parameters but you provided %f"%len(parlist)
    elif len(parlist)<6:
        print "to few parameters, setting p1 to p3 in the exponential cutoff to 1"
        for i in range(6-len(parlist)):
            parlist.append(1)
    if x< par[3]:
        return parlist[0]*(x/parlist[1])**parlist[2]
    else:
        return parlist[0]*(x/parlist[1])**parlist[2]*np.exp(-((x-parlist[3])/parlist[4]+par[5]*np.log(x/parlistp[3])+par[6]*(log(x/parlist[3]))**2))
   '''     
def test(x,parlist=[1,2]):
    if len(parlist)!=2:
        print "function test needs 2 arguments but %f given"%len(parlist)
    return 10**(parlist[0]+parlist[1]*x)
    
