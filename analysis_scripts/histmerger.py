import ROOT
import numpy as np
import sys
import os
import random
import SndlhcGeo
from tqdm import tqdm
from statistics import mode
from argparse import ArgumentParser
parser=ArgumentParser()
parser.add_argument("--l",dest="legend",nargs="*")
parser.add_argument("--h",dest="histname")
parser.add_argument("--o",dest="outfile",default="outfile")
parser.add_argument("--f",dest="histfiles",nargs="*")
parser.add_argument("--allhists",action='store_true',dest='allhists')
parser.add_argument("--t",dest="title",default="title")
args=parser.parse_args()
hstacks={}
c=ROOT.TCanvas()
files={}

hists={}
if args.allhists:
	hstack=ROOT.THStack(args.histfiles[0],"stack of histograms from {}".format(args.histfiles[0]))
	file=ROOT.TFile(args.histfiles[0],'read')
	for i,key in enumerate(file.GetListOfKeys()):
		hists[i]=file.Get(key.GetName())		
		hists[i].Scale(1./hists[i].Integral())
		hists[i].SetLineColor(2+i)
		hists[i].SetLineStyle(1)
		hists[i].SetOption("B")
		hists[i].SetTitle(key.GetName())
		hists[i].SetFillColorAlpha(2+i,0.4)
		hists[i].SetFillStyle(1001)
		hstack.Add(hists[i].Clone("qdcus{}".format(i)))
else:
	hstack=ROOT.THStack(args.histname,"stack of {} histograms".format(args.histname))
	for i,filename in enumerate(args.histfiles):	
		#each dict entry is a stack 
		files[i]=ROOT.TFile(filename,'read')
		hists[i]=files[i].Get(args.histname)		
		hists[i].Scale(1./hists[i].Integral())
		hists[i].SetLineColor(2+i)
		hists[i].SetLineStyle(1)
		hists[i].SetOption("B")
		if len(args.legend)>i:
			hists[i].SetTitle(args.legend[i])
		else:
			hists[i].SetTitle(filename.replace(".root",""))
		hists[i].SetFillColorAlpha(2+i,0.4)
		hists[i].SetFillStyle(1001)
		hstack.Add(hists[i].Clone("qdcus{}".format(i)))
#print(hstack.GetHists().GetSize())
hstack.Draw("HIST nostack")
hstack.GetXaxis().SetTitle(hists[0].GetXaxis().GetTitle())
hstack.GetYaxis().SetTitle(hists[0].GetYaxis().GetTitle())
c.Update()
c.BuildLegend(0.7,0.75,0.95,0.95)
c.SaveAs("{}.png".format(args.outfile))			

'''for histname in [key.GetName() for key in file.GetListOfKeys()]:
			if histname in hstacks.keys():
				hstacks[histname].Add(file.Get(histname))
			else:
				hstacks[histname]=ROOT.THStack(histname,"stack of {} histograms".format(histname))
				hstacks[histname].Add(file.Get(histname))'''

'''with ROOT.TFile('{}.root'.format(args.outfile),'RECREATE') as ofile:
	for stack in hstacks:
		ofile.WriteObject(hstacks[stack],stack)'''