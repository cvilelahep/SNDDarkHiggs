import ROOT
import numpy as np
import math
import re
from tqdm import tqdm
from argparse import ArgumentParser
parser=ArgumentParser()
parser.add_argument("--o",dest="outfile",default="outfile")
parser.add_argument("--f",dest="histfiles",nargs="*")
args=parser.parse_args()
hstacks={}
c=ROOT.TCanvas()
hists={}
fitfunc={}
mu={}
ids={}
#def cbf(x,p):
#    A=((p[2]/p[1])**p[2])*math.exp(-p[1]**2/2)
#    B=p[2]/p[1]-p[1]
#    C=p[2]/p[1]*1/(p[2]-1)*math.exp(-p[1]**2/2)
#    D=math.sqrt(math.pi/2)*(1+math.erf(p[1]/math.sqrt(2)))
#    N=1/(C+D)
#    if x-p[0]>p[1]:
#        return N*A*(B+(x-p[0]))**(-p[2])
#    if x-p[0]<=p[1]:
#        return N*math.exp(-(x-p[0])**2/2)
#print(cbf(25,[25,1,2]))
#rootcbf=ROOT.TF1("cbf",cbf,1,250,3)
#rootcbf.SetParLimits(1,1,250)
#rootcbf.SetParLimits(2,2,5)
f = open(f"{args.outfile}.txt", "w")
f.write("Now the file has more content!")

file=ROOT.TFile(args.histfiles[0],'read')
for i,key in enumerate(file.GetListOfKeys()):
    hists[i]=file.Get(key.GetName())    	
    mu[i]=hists[i].GetMaximumBin()*2
    #hists[i].Fit("gaus")
    #histgaus[i]=hists[i].GetFunction("gaus")
    #rootcbf.SetParameters(hists[i].GetMean(),1,2)
    #mu[i]=fitfunc[i].GetParameter(1)
    
    id=re.search(r'plane(?P<plane>\d+).*?bar(?P<bar>\d+)',hists[i].GetTitle())
    id.groupdict()    
    ids[i]=int(id["plane"])*10+int(id["bar"])
avgmu=0
for key in mu.keys():   
    avgmu+=mu[key]/len(mu.keys())
f.write(f'average mu to align to: {avgmu}')
for i in range(len(ids)):
    f.write("shiftdict[{}]={}".format(ids[i],mu[i]-avgmu))
f.close()
