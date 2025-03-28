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
parser.add_argument("--g",dest="graphname")
parser.add_argument("--o",dest="outfile",default="outfile")
parser.add_argument("--f",dest="graphfiles",nargs="*",required=True)
parser.add_argument("--allgraphs",action='store_true',dest='allgraphs')
parser.add_argument("--log",action='store_true',dest='logscale')
args=parser.parse_args()
hstacks={}
c=ROOT.TCanvas()
files={}
graphs={}
if args.allgraphs:
	gstack=ROOT.TMultiGraph(args.graphfiles[0],"stack of graphs from {}".format(args.graphfiles[0]))
	file=ROOT.TFile(args.graphfiles[0],'read')
	for i,key in enumerate(file.GetListOfKeys()):
		graphs[i]=file.Get(key.GetName())	
		graphs[i].SetLineColor(2+i)
		graphs[i].SetLineStyle(1)
		graphs[i].Draw("AP")
		if args.legend:
			if len(args.legend)>i:
				graphs[i].SetTitle(args.legend[i])
		else:
			graphs[i].SetTitle(args.graphfiles[0].replace(".root",""))
		gstack.Add(graphs[i])	
else:
	gstack=ROOT.TMultiGraph(args.graphname,"stack of {} graphs".format(args.graphname))
	gstack.SetTitle("{}".format(args.graphname))
	for i,filename in enumerate(args.graphfiles):	
		#each dict entry is a stack 
		files[i]=ROOT.TFile(filename,'read')
		graphs[i]=files[i].Get(args.graphname)	
		graphs[i].SetLineColor(2+i)
		graphs[i].SetLineStyle(1)
		graphs[i].Draw("AP")
		if args.legend:
			if len(args.legend)>i:
				graphs[i].SetTitle(args.legend[i])
		else:
			graphs[i].SetTitle(filename.replace(".root",""))
		gstack.Add(graphs[i])
if args.logscale:
	c.SetLogx()
#print(hstack.Getgraphs().GetSize())
gstack.Draw("ALP")
c.Update()
c.BuildLegend(0.7,0.75,0.95,0.95)
with ROOT.TFile('{}.root'.format(args.outfile),'RECREATE') as ofile:
	ofile.WriteObject(gstack,"merged_{}".format(args.outfile))
c.SaveAs("{}.png".format(args.outfile))			

'''for graphname in [key.GetName() for key in file.GetListOfKeys()]:
			if graphname in hstacks.keys():
				hstacks[graphname].Add(file.Get(graphname))
			else:
				hstacks[graphname]=ROOT.THStack(graphname,"stack of {} histograms".format(graphname))
				hstacks[graphname].Add(file.Get(graphname))'''

'''with ROOT.TFile('{}.root'.format(args.outfile),'RECREATE') as ofile:
	for stack in hstacks:
		ofile.WriteObject(hstacks[stack],stack)'''