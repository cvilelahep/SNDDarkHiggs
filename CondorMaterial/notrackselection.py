import ROOT
import numpy as np
import sys
import os
import random
import SndlhcGeo
from tqdm import tqdm
from statistics import mode
from argparse import ArgumentParser
from array import array
parser=ArgumentParser()
parser.add_argument("--realdata",action='store_true',dest='realdata')
parser.add_argument("--o",dest="outfile",default="outfile")
parser.add_argument("--f",dest="digfile",nargs="*")
parser.add_argument("--dir",dest="rootdir",default=".")
parser.add_argument("--usingdir",action='store_true',dest='usingdir')

args=parser.parse_args()


addlist=[]
if args.realdata:
	treeinfoname="rawConv"
else:
	treeinfoname="cbmsim"
if args.usingdir:
	treedigi=ROOT.TChain(treeinfoname)
	for subdir in [subdir for subdir in os.listdir(args.rootdir) if os.path.isdir("{}/{}".format(args.rootdir,subdir))]:
		for file in [filename for filename in os.listdir('{}/{}'.format(args.rootdir,subdir)) if filename.endswith("_dig.root")]:
			treedigi.Add("{}/{}/{}".format(args.rootdir,subdir,file))
elif len(args.digfile)>=1:
	treedigi=ROOT.TChain(treeinfoname)
	for file in args.digfile:
		treedigi.Add(file)	
if not treedigi:
	print("failed to get treedigi, exiting")
	sys.exit()


geoname=[filename for filename in os.listdir('.') if filename.startswith("geofile")]
geo=SndlhcGeo.GeoInterface(geoname[0])	
scifi=geo.modules['Scifi']		

prevtime=0
for i_event, eventdigi in tqdm(enumerate(treedigi)) :
	reject=False
	realdataaproved=True
	veto=False
	if args.realdata:
		if eventdigi.EventHeader.GetBeamMode()!=11 or eventdigi.EventHeader.isIP1()!=True or eventdigi.EventHeader.GetEventTime()-prevtime<100:
			realdataaproved=False
	prevtime=eventdigi.EventHeader.GetEventTime()
	for hit in eventdigi.Digi_MuFilterHits:
		if hit.GetSystem()==1:
			veto=True
			break

#flag decides if veto cut is applied (no veto hits)
	if realdataaproved:
		tempscifiplane=[]	
		tempmufiupplane=[]
		scifipassx,scifipassy=False,False
		for hit in eventdigi.Digi_ScifiHits:
			if not hit.isValid():
				continue
						

			#print('a:{},{},{}'.format(a.X(),a.Y(),a.Z()))
			#print('b:{},{},{}'.format(b.X(),b.Y(),b.Z()))							
			#print('track:{},{},{}\n\n'.format(getTrackAtZ(track,a.Z()).X(),getTrackAtZ(track,a.Z()).Y(),getTrackAtZ(track,a.Z()).Z()))
			if hit.isVertical():
				if hit.GetStation()==5:
					scifipassx=True							
			else:
				if hit.GetStation()==5:
					scifipassy=True
		if not(scifipassy and scifipassx):
			reject=True
			continue
		for hit in eventdigi.Digi_MuFilterHits:
			if not hit.isValid():
				continue
			if hit.GetSystem()==2:
					tempmufiupplane.append(hit.GetPlane())
		if tempmufiupplane:
			if len(tempmufiupplane)>5 or tempmufiupplane.count(mode(tempmufiupplane))>2:
				reject=True
				continue
		elif not tempmufiupplane:
			reject=True
			continue
		if not reject:
			addlist.append(i_event)
	#if len(addlist)>3:
	#	break
print("adding {} events to new file".format(len(addlist)))
filedigi=treedigi.GetFile()
with ROOT.TFile("{}.root".format(args.outfile),'RECREATE') as newfiledigi:
	for key in filedigi.GetListOfKeys():
		#if key.GetName()!=treeinfoname:
		if key.GetName()!=treeinfoname and key.GetName()!="cbmsim":
			newfiledigi.WriteObject(filedigi.Get(key.GetName()),key.GetName())
	newtreedigi=treedigi.CloneTree(0)
	for i in addlist:
		treedigi.GetEntry(i)
		newtreedigi.Fill()
	#newfiledigi.WriteObject(newtreedigi,treeinfoname)
	newfiledigi.WriteObject(newtreedigi,treeinfoname)
