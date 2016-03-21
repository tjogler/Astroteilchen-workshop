import numpy as np

def func(name):
    if name=="pl":
        return [pl,3,[1,0,1]]
    elif name=="bpl":
        return [bpl,4,[1,0,0,1]]
    elif name=="plcut":
        return [plcut,4,[1,0,1,1]]
    elif name=="test":
        return [test,2,[1,1]]


def pl(x,parlist):
    

def bpl(x,parlist):
    

def plcut(x,par):
    

def test(x,parlist=[1,2]):
    if len(parlist)!=2:
        print "function test needs 2 arguments but %f given"%len(parlist)
    return 10**(parlist[0]+parlist[1]*x)
    
