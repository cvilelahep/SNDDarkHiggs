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
parser.add_argument("--SciFicut",action='store_true',dest='scificut')
parser.add_argument("--UScut",action='store_true',dest='uscut')
parser.add_argument("--Vetocut",action='store_true',dest='vetocut')
parser.add_argument("--Anglecut",action='store_true',dest='anglecut')
parser.add_argument("--XYcut",action='store_true',dest='xycut')
parser.add_argument("--l",dest="legend",default="placeholder")
parser.add_argument("--o",dest="outfile",default="outfile")
parser.add_argument("--f",dest="digfile",nargs="*")
parser.add_argument("--pfunc",dest="pfile",required=True)
parser.add_argument("--dir",dest="rootdir",default=".")
parser.add_argument("--usingdir",action='store_true',dest='usingdir')
parser.add_argument("--w",action='store_true',dest='weighted')
parser.add_argument("--MC",action='store_true',dest='MCtest')
parser.add_argument("--shift",action='store_true',dest='shifted')
args=parser.parse_args()

shiftdict={}
shiftdict[4]=-16.857142857142875
shiftdict[14]=11.142857142857125
shiftdict[24]=17.142857142857125
shiftdict[34]=-22.857142857142875
shiftdict[44]=31.142857142857125
shiftdict[6]=-18.857142857142875
shiftdict[16]=-10.857142857142875
shiftdict[26]=17.142857142857125
shiftdict[36]=-26.857142857142875
shiftdict[46]=27.142857142857125
shiftdict[2]=7.142857142857125
shiftdict[12]=-4.857142857142875
shiftdict[22]=11.142857142857125
shiftdict[32]=-12.857142857142875
shiftdict[42]=35.142857142857125
shiftdict[5]=-12.857142857142875
shiftdict[15]=15.142857142857125
shiftdict[25]=15.142857142857125
shiftdict[35]=-8.857142857142875
shiftdict[45]=15.142857142857125
shiftdict[13]=7.142857142857125
shiftdict[23]=1.1428571428571246
shiftdict[33]=-18.857142857142875
shiftdict[43]=13.142857142857125
shiftdict[1]=-10.857142857142875
shiftdict[3]=-12.857142857142875
shiftdict[7]=9.142857142857125
shiftdict[17]=-6.857142857142875
shiftdict[27]=9.142857142857125
shiftdict[37]=-26.857142857142875
shiftdict[47]=-4.857142857142875
shiftdict[11]=-8.857142857142875
shiftdict[21]=-2.8571428571428754
shiftdict[31]=-24.857142857142875
shiftdict[41]=11.142857142857125

pfile=ROOT.TFile(args.pfile,'read')
p1=pfile.Get("p1")
p2=pfile.Get("p2")
if args.realdata:
	treeinfoname="rawConv"
else:
	treeinfoname="cbmsim"
if args.usingdir:
	treedigi=ROOT.TChain(treeinfoname)
	treereco=ROOT.TChain(treeinfoname)
	for subdir in [subdir for subdir in os.listdir(args.rootdir) if os.path.isdir("{}/{}".format(args.rootdir,subdir))]:
		for file in [filename for filename in os.listdir('{}/{}'.format(args.rootdir,subdir)) if filename.endswith("_dig.root")]:
			treedigi.Add("{}/{}/{}".format(args.rootdir,subdir,file))
			treereco.Add("{}/{}/{}".format(args.rootdir,subdir,file).replace('.root','__muonReco.root'))
elif len(args.digfile)>=1:
	treedigi=ROOT.TChain(treeinfoname)
	treereco=ROOT.TChain(treeinfoname)
	for file in args.digfile:
		treedigi.Add(file)
		treereco.Add(file.replace('.root','__muonReco.root'))		
if not treedigi:
	print("failed to get treedigi, exiting")
	sys.exit()
if not treereco:
	print("failed to get treereco, exiting")
	sys.exit()
eventcount=0
passcount=0
passvetocount=0
wpasscount=0
weventcount=0
xarray=[]
yarray=[]
angleXZ=[]
angleYZ=[]
minqdcdw=[]
muonhits=[]
muonhitsbar=[]
eventid=[]
hittrackmaxdiff=[]
totalscifihits=[]
maxscifihits=[]
totalmufiuphits=[]
maxmufiuphits=[]
totalmufidwhits=[]
maxmufidwhits=[]
bestqdcdistup=[]
hittheveto=[]
noUS=0
vetoinit=[0,0]
vetoxycut=[0,0]
vetoUScut=[0,0]
vetoscificut=[0,0]
vetoanglecut=[0,0]
def getTrackAtZ(ftrack,z = 280):
	start=ftrack.getStart()
	mom=ftrack.getTrackMom()
	targetZ=z
	#get k from the z coordinates
	k=(start.Z()-targetZ)/mom.Z()
	#calculate X and Y and return TVector3
	targetX=start.X()-k*mom.X()
	targetY=start.Y()-k*mom.Y()
	target= ROOT.TVector3(targetX,targetY,targetZ)
	return target
geoname=[filename for filename in os.listdir('.') if filename.startswith("geofile")]
geo=SndlhcGeo.GeoInterface(geoname[0])	
scifi=geo.modules['Scifi']	
mufilter=geo.modules['MuFilter']	

prevtime=0
for i_event, [eventdigi,eventreco] in tqdm(enumerate(zip(treedigi,treereco))) :
	reject=False
	realdataaproved=True
	a=ROOT.TVector3()
	b=ROOT.TVector3()
	veto=False
	eventcount+=1
	if args.realdata:
		if eventdigi.EventHeader.GetBeamMode()!=11 or eventdigi.EventHeader.isIP1()!=True or eventdigi.EventHeader.GetEventTime()-prevtime<100:
			realdataaproved=False
	prevtime=eventdigi.EventHeader.GetEventTime()
	for hit in eventdigi.Digi_MuFilterHits:
		if hit.GetSystem()==1:
			veto=True
			break
	if args.MCtest and random.random()<0.2:
		veto=False 			
#flag decides if veto cut is applied (no veto hits)
	if realdataaproved:
		tempscifiplane=[]
		if eventreco.Reco_MuonTracks.GetEntries()>0:
			for track in eventreco.Reco_MuonTracks:
				tempqdcmindw=[]
				tempqdcminup=[]
				tempmufiupplane=[]
				bestqdcperplaneup={} 
				minhittrackdiffperplaneup={}
				bestqdcperplanedw={} 
				minhittrackdiffperplanedw={}
				IDofminhitup={}
				IDofminhitdw={}
				distofbestqdcup={}
				associatedbar={}
				hitcountscifi=0
				maxhittrackdiffx=0
				maxhittrackdiffy=0
				if veto and (not reject):
					vetoinit[1]+=1
				if (not veto) and (not reject):
					vetoinit[0]+=1
				
				vetocoords=getTrackAtZ(track)
				if args.xycut and not (np.logical_and(-44<vetocoords.X()<-10,18<vetocoords.Y()<52.5) and np.logical_and(-47<getTrackAtZ(track,500).X()<-8,15.5<getTrackAtZ(track,500).Y()<54.5)): 
					reject=True
					continue
				if veto and (not reject):
					vetoxycut[1]+=1
				if (not veto) and (not reject):
					vetoxycut[0]+=1
				if args.scificut:
					scifipassx,scifipassy=False,False		
					for hit in eventdigi.Digi_ScifiHits:
						if not hit.isValid():
							continue
							
						scifi.GetSiPMPosition(hit.GetChannelID(),a,b)
						#print('a:{},{},{}'.format(a.X(),a.Y(),a.Z()))
						#print('b:{},{},{}'.format(b.X(),b.Y(),b.Z()))							
						#print('track:{},{},{}\n\n'.format(getTrackAtZ(track,a.Z()).X(),getTrackAtZ(track,a.Z()).Y(),getTrackAtZ(track,a.Z()).Z()))
						if hit.isVertical():
							if hit.GetStation()==5:
								scifipassx=True
							if maxhittrackdiffx<(getTrackAtZ(track,a.Z()).X()-(a.X()/2+b.X()/2))**2:
								maxhittrackdiffx=(getTrackAtZ(track,a.Z()).X()-(a.X()/2+b.X()/2))**2
							
						else:
							if hit.GetStation()==5:
								scifipassy=True
							if maxhittrackdiffy<(getTrackAtZ(track,a.Z()).Y()-(a.Y()/2+b.Y()/2))**2:
								maxhittrackdiffy=(getTrackAtZ(track,a.Z()).Y()-(a.Y()/2+b.Y()/2))**2
					if not(scifipassy and scifipassx and max(maxhittrackdiffx,maxhittrackdiffy)<5):
						reject=True
						continue
				if veto and (not reject):
					vetoscificut[1]+=1
				if (not veto) and (not reject):
					vetoscificut[0]+=1
				for hit in eventdigi.Digi_MuFilterHits:
					if not hit.isValid():
						continue
					if hit.GetSystem()==1:
							continue
					elif hit.GetSystem()==2:
							tempmufiupplane.append(hit.GetPlane())
							vecA=ROOT.TVector3()
							vecB=ROOT.TVector3()
							mufilter.GetPosition(hit.GetDetectorID(),vecA,vecB)
							if hit.isVertical():
								'''if abs(getTrackAtZ(track,vecA.Z()).X()-(vecA.X()/2+vecB.X()/2))<6:	
									totalsignal=0
									for [sipm,qdc] in hit.GetAllSignals():
										if not hit.isShort(sipm):
											totalsignal+=qdc
									if hit.GetPlane() not in bestqdcperplaneup.keys():
										bestqdcperplaneup[hit.GetPlane()]=totalsignal
										IDofminhitup[hit.GetPlane()]=hit.GetDetectorID()
										distofbestqdcup[hit.GetPlane()]=abs(getTrackAtZ(track,vecA.Z()).X()-(vecA.X()/2+vecB.X()/2))
									elif totalsignal>bestqdcperplaneup[hit.GetPlane()]:
										bestqdcperplaneup[hit.GetPlane()]=totalsignal
										IDofminhitup[hit.GetPlane()]=hit.GetDetectorID()
										distofbestqdcup[hit.GetPlane()]=abs(getTrackAtZ(track,vecA.Z()).X()-(vecA.X()/2+vecB.X()/2))'''
								if hit.GetPlane() not in minhittrackdiffperplaneup or (minhittrackdiffperplaneup[hit.GetPlane()]>(getTrackAtZ(track,vecA.Z()).X()-(vecA.X()/2+vecB.X()/2))**2):
									minhittrackdiffperplaneup[hit.GetPlane()]=(getTrackAtZ(track,vecA.Z()).X()-(vecA.X()/2+vecB.X()/2))**2
									totalsignal=0
									for [sipm,qdc] in hit.GetAllSignals():
										if not hit.isShort(sipm):
											totalsignal+=qdc
									bestqdcperplaneup[hit.GetPlane()]=totalsignal
									IDofminhitup[hit.GetPlane()]=hit.GetDetectorID()
									distofbestqdcup[hit.GetPlane()]=abs(getTrackAtZ(track,vecA.Z()).X()-(vecA.X()/2+vecB.X()/2))
									associatedbar[hit.GetPlane()]=hit.GetDetectorID()%1000
							else:
								'''if abs(getTrackAtZ(track,vecA.Z()).Y()-(vecA.Y()/2+vecB.Y()/2))<6:
									totalsignal=0
									for [sipm,qdc] in hit.GetAllSignals():
										if not hit.isShort(sipm):
											totalsignal+=qdc
									if hit.GetPlane() not in bestqdcperplaneup.keys():
										bestqdcperplaneup[hit.GetPlane()]=totalsignal
										IDofminhitup[hit.GetPlane()]=hit.GetDetectorID()
										distofbestqdcup[hit.GetPlane()]=abs(getTrackAtZ(track,vecA.Z()).Y()-(vecA.Y()/2+vecB.Y()/2))
									elif totalsignal>bestqdcperplaneup[hit.GetPlane()]:
										bestqdcperplaneup[hit.GetPlane()]=totalsignal
										IDofminhitup[hit.GetPlane()]=hit.GetDetectorID()
										distofbestqdcup[hit.GetPlane()]=abs(getTrackAtZ(track,vecA.Z()).Y()-(vecA.Y()/2+vecB.Y()/2))'''
								if hit.GetPlane() not in minhittrackdiffperplaneup or (minhittrackdiffperplaneup[hit.GetPlane()]>(getTrackAtZ(track,vecA.Z()).Y()-(vecA.Y()/2+vecB.Y()/2))**2):
									minhittrackdiffperplaneup[hit.GetPlane()]=(getTrackAtZ(track,vecA.Z()).Y()-(vecA.Y()/2+vecB.Y()/2))**2
									totalsignal=0
									for [sipm,qdc] in hit.GetAllSignals():
										if not hit.isShort(sipm):
											totalsignal+=qdc
									bestqdcperplaneup[hit.GetPlane()]=totalsignal
									IDofminhitup[hit.GetPlane()]=hit.GetDetectorID()
									distofbestqdcup[hit.GetPlane()]=abs(getTrackAtZ(track,vecA.Z()).Y()-(vecA.Y()/2+vecB.Y()/2))
									associatedbar[hit.GetPlane()]=hit.GetDetectorID()%1000
					elif hit.GetSystem()==3:
							continue
							vecA=ROOT.TVector3()
							vecB=ROOT.TVector3()
							hit.GetPosition(vecA,vecB)
							if hit.isVertical():
								if abs(getTrackAtZ(track,vecA.Z()).X()-(vecA.X()/2+vecB.X()/2))<6:
									totalsignal=0
									for [sipm,qdc] in hit.GetAllSignals():
										if not hit.isShort(sipm):
											totalsignal+=qdc
									if hit.GetPlane() not in bestqdcperplanedw.keys():
										bestqdcperplanedw[hit.GetPlane()]=totalsignal
										IDofminhitdw[hit.GetPlane()]=hit.GetDetectorID()
									elif totalsignal>bestqdcperplanedw[hit.GetPlane()]:
										bestqdcperplanedw[hit.GetPlane()]=totalsignal
										IDofminhitdw[hit.GetPlane()]=hit.GetDetectorID()
								'''if hit.GetPlane() not in minhittrackdiffperplanedw or (minhittrackdiffperplanedw[hit.GetPlane()]>(getTrackAtZ(track,vecA.Z()).X()-(vecA.X()/2+vecB.X()/2))**2):
									minhittrackdiffperplanedw[hit.GetPlane()]=(getTrackAtZ(track,vecA.Z()).X()-(vecA.X()/2+vecB.X()/2))**2
									totalsignal=0
									for [sipm,qdc] in hit.GetAllSignals():
										if not hit.isShort(sipm):
											totalsignal+=qdc
									bestqdcperplanedw[hit.GetPlane()]=totalsignal'''
									
							else:
								if abs(getTrackAtZ(track,vecA.Z()).Y()-(vecA.Y()/2+vecB.X()/2))<6 :
									totalsignal=0
									for [sipm,qdc] in hit.GetAllSignals():
										if not hit.isShort(sipm):
											totalsignal+=qdc
									if hit.GetPlane() not in bestqdcperplanedw:
										bestqdcperplanedw[hit.GetPlane()]=totalsignal
										IDofminhitdw[hit.GetPlane()]=hit.GetDetectorID()
									elif totalsignal>bestqdcperplanedw[hit.GetPlane()]:
										bestqdcperplanedw[hit.GetPlane()]=totalsignal
										IDofminhitdw[hit.GetPlane()]=hit.GetDetectorID()
								'''if hit.GetPlane() not in minhittrackdiffperplanedw or (minhittrackdiffperplanedw[hit.GetPlane()]>(getTrackAtZ(track,vecA.Z()).Y()-(vecA.Y()/2+vecB.Y()/2))**2):
									minhittrackdiffperplanedw[hit.GetPlane()]=(getTrackAtZ(track,vecA.Z()).Y()-(vecA.Y()/2+vecB.Y()/2))**2
									totalsignal=0
									for [sipm,qdc] in hit.GetAllSignals():
										if not hit.isShort(sipm):
											totalsignal+=qdc
									bestqdcperplanedw[hit.GetPlane()]=totalsignal'''
					
					else:
						sys.exit("Hit must have a system")
				if args.uscut:
					if tempmufiupplane:
						if len(tempmufiupplane)>5 or tempmufiupplane.count(mode(tempmufiupplane))>2:
							reject=True
							continue
					elif not tempmufiupplane:
						reject=True
						continue
				if veto and (not reject):
					vetoUScut[1]+=1
				if (not veto) and (not reject):
					vetoUScut[0]+=1
				if args.anglecut:
					if abs(track.getAngleYZ())>0.02 or abs(track.getAngleXZ())>0.01:
						reject=True
						continue
				if veto and (not reject):
					vetoanglecut[1]+=1
				if (not veto) and (not reject):
					vetoanglecut[0]+=1
				if bestqdcperplaneup and (not reject):
					muonhits.append(bestqdcperplaneup)
					muonhitsbar.append(associatedbar)
					if veto and (not reject):
						passvetocount+=1
						hittheveto.append(1)
					if (not veto) and (not reject):
						hittheveto.append(0)
	if len(muonhits)!=len(hittheveto):
		print("Error: muonhits and hittheveto have different lengths")
		exit()			

	#if len(muonhits)>10:
	#	break			

#A=No veto likelihood<2
#B=Veto likelihood<2
#C=No veto likelihood>2
#D=Veto likelihood>2
print("vetoinit: {} hit, {} nohit".format(vetoinit[1],vetoinit[0]))
print("vetoxycut: {} hit, {} nohit, cut was {}".format(vetoxycut[1],vetoxycut[0],args.xycut))
print("vetoscificut: {} hit, {} nohit, cut was {}".format(vetoscificut[1],vetoscificut[0],args.scificut))
print("vetoUScut: {} hit, {} nohit, cut was {}".format(vetoUScut[1],vetoUScut[0],args.uscut))
print("vetoanglecut: {} hit, {} nohit, cut was {}".format(vetoanglecut[1],vetoanglecut[0],args.anglecut))
f=open("ABCD_{}.txt".format(args.outfile),'w')
f.write("vetoinit: {} hit, {} nohit\n".format(vetoinit[1],vetoinit[0]))
f.write("vetoxycut: {} hit, {} nohit, cut was {}\n".format(vetoxycut[1],vetoxycut[0],args.xycut))
f.write("vetoscificut: {} hit, {} nohit, cut was {}\n".format(vetoscificut[1],vetoscificut[0],args.scificut))
f.write("vetoUScut: {} hit, {} nohit, cut was {}\n".format(vetoUScut[1],vetoUScut[0],args.uscut))
f.write("vetoanglecut: {} hit, {} nohit, cut was {}\n".format(vetoanglecut[1],vetoanglecut[0],args.anglecut))

Abkg_list=array('d')
Abkg_error=array('d')
Azone_list=array('d')
Azone_error=array('d')
threshold_list=array('d')
logcountlist=[]
f.write("\n no-veto approved \n")
for i in range(len(muonhits)):
	logcount=0
	for plane in muonhits[i].keys():
		muonhit_qdc=muonhits[i][plane]-shiftdict[int(muonhitsbar[i][plane])+int(plane)*10] if (args.shifted and (muonhitsbar[i][plane])+int(plane)*10 in shiftdict.keys()) else muonhits[i][plane]
		if args.shifted and not ((muonhitsbar[i][plane])+int(plane)*10 in shiftdict.keys()):
			continue
		logcount+=np.log(p1.GetBinContent(int(int(muonhit_qdc)/2)+1)) if p1.GetBinContent(int(int(muonhit_qdc)/2)+1)!=0 else 0
		logcount-=np.log(p2.GetBinContent(int(int(muonhit_qdc)/2)+1)) if p2.GetBinContent(int(int(muonhit_qdc)/2)+1)!=0 else 0
	logcountlist.append(logcount)
	if hittheveto[i]==0:
		f.write("logcount={}\n".format(logcount))
for logcount_threshold in tqdm(np.arange(-10,10,0.5)):
	B,C,D,A_zone=0.,0.,0.,0.
	threshold_list.append(logcount_threshold)
	D=len([logcount for logcount,veto in zip(logcountlist,hittheveto) if (logcount>=logcount_threshold and veto==1)])
	B=len([logcount for logcount,veto in zip(logcountlist,hittheveto) if (logcount<logcount_threshold and veto==1)])
	C=len([logcount for logcount,veto in zip(logcountlist,hittheveto) if (logcount>=logcount_threshold and veto==0)])
	A_zone=len([logcount for logcount,veto in zip(logcountlist,hittheveto) if (logcount<logcount_threshold and veto==0)])
	Abkg=(C/D)*B if D!=0 else -1
	Abkg_list.append(Abkg)
	Abkg_error.append(np.sqrt(Abkg))
	Azone_list.append(A_zone)
	Azone_error.append(0)
	f.write("\n threshold={}: \n".format(logcount_threshold))
	f.write("from approved {} events, {} hit veto\n".format(len(muonhits),len([i for i in hittheveto if i==1])))
	f.write('regions with only background: B (veto log<{})={},C(no veto log>{})={},D(veto log>{})={}\n'.format(logcount_threshold,B,logcount_threshold,C,logcount_threshold,D))
	f.write('expected background in signal region: A={}\n , total {} signal region events\n'.format(Abkg,A_zone))
#hAzone=ROOT.TGraph(len(threshold_list),threshold_list,Azone_list)
#hAbkg=ROOT.TGraphErrors(len(threshold_list),threshold_list,Abkg_list,Azone_error,Abkg_error)
#with ROOT.TFile("ABCD_{}.root".format(args.outfile),'RECREATE') as ofile:
	#ofile.WriteObject(hAbkg,'Abkg')
	#ofile.WriteObject(hAzone,'Azone')
f.close()


	
print('process ended successfully')
