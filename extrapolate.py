#!/usr/bin/env python

import os,sys,glob
import pyfits
import numpy as np
import ROOT as root
from math import degrees
import argparse

def readFile(infile):
    
    fin=open(infile)
    lines=fin.readlines()
    
    graph=root.TGraphErrors()
    print lines

    for l in lines:
        v=l.split(' ')
        graph.SetPoint(lines.index(l),float(v[0]),float(v[1]))
        graph.SetPointError(lines.index(l),0,float(v[2]))

    return graph

def extrapolate(infile):
    gr=readFile(infile)
    can=root.TCanvas("graph","graph",800,600)
    can.cd()
    gr.Draw("PA")
    var=raw_input("Press return to exit")

parser=argparse.ArgumentParser()
parser.add_argument("--file",help="file with the data")
args=parser.parse_args()
extrapolate(args.file)
