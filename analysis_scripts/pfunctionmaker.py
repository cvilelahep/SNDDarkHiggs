import ROOT
import numpy as np
import sys
import os
import SndlhcGeo
from tqdm import tqdm
from statistics import mode
from argparse import ArgumentParser
parser=ArgumentParser()
parser.add_argument("--realdata",action='store_true',dest='realdata')
parser.add_argument("--l",dest="legend",default="placeholder")
parser.add_argument("--o",dest="outfile",default="outfile")
parser.add_argument("--f",dest="digfile",nargs="*")
parser.add_argument("--dir",dest="rootdir",default=".")
parser.add_argument("--usingdir",action='store_true',dest='usingdir')
parser.add_argument("--w",action='store_true',dest='weighted')
parser.add_argument("--shift",action='store_true',dest='shifted')
args=parser.parse_args()

#This is a shift dictionary, that offsets the QDC of each US bar to the average most probable QDC
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


#Importing files
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
noUS=0
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
				init+=1
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
				#calculate angle for anglecut selection efficiency
				station_stats = {}  # To store info about visited stations
				for hit in eventdigi.Digi_ScifiHits:						
					if not hit.isValid():
						continue							
					scifi.GetSiPMPosition(hit.GetChannelID(),a,b)
					# Track station stats
					station=hit.GetStation()
					if station not in station_stats:
						station_stats[station] = {
							"sum_x": 0,
							"sum_y": 0,
							"sum_z": 0,	
							"countx": 0,
							"county": 0,
							"countz": 0,
							"has_vertical": False,
							"has_horizontal": False,
						}
					if hit.isVertical():
						station_stats[station]["sum_x"]+=a.X()/2+b.X()/2
						station_stats[station]["countx"] += 1
						station_stats[station]["sum_z"]+=a.Z()/2+b.Z()/2
						station_stats[station]["countz"] += 1
						station_stats[station]["has_vertical"] = True
					else:
						station_stats[station]["sum_y"]+=a.Y()/2+b.Y()/2
						station_stats[station]["county"] += 1
						station_stats[station]["sum_z"]+=a.Z()/2+b.Z()/2
						station_stats[station]["countz"] += 1
						station_stats[station]["has_horizontal"] = True							
				# Find the lowest station with both vertical and horizontal hits
				bothvandh=False#flag that signals if at least one station has both vertical and horizontal hits
				for station in sorted(station_stats.keys()):
					stats = station_stats[station]
					if stats["has_vertical"] and stats["has_horizontal"]:
						bothvandh=True
						avg_x = stats["sum_x"] / stats["countx"]
						avg_y = stats["sum_y"] / stats["county"]
						avg_z = stats["sum_z"] / stats["countz"]
						break
				trackangleXZ= (getTrackAtZ(track,500).X()-avg_x)/ (500-avg_z) if bothvandh else track.getAngleXZ()
				trackangleYZ= (getTrackAtZ(track,500).Y()-avg_y)/ (500-avg_z) if bothvandh else track.getAngleYZ()	
				
				vetocoords=ROOT.TVector3(getTrackAtZ(track,500).X()-(500-280)*trackangleXZ,getTrackAtZ(track,500).Y()-(500-280)*trackangleYZ,280) #280 is the veto Z
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
					else:
						sys.exit("Hit must have a system")
				if bestqdcperplaneup:
					muonhits.append(bestqdcperplaneup)
					muonhitsbar.append(associatedbar)
	#if len(muonhits)>20:
	#	break
#hqdcmax=ROOT.TH1D("max_qdc","Maximum QDC for each event for {} events;maxqdc;count".format(args.legend),200,0,200)
p1func=ROOT.TH1D("p1func","p1func;hit_qdc;\%",250,0,500)
p2func=ROOT.TH1D("p2func","p2func;hit_qdc;\%",250,0,500)
for i in range(len(muonhits)):
	for plane in muonhits[i].keys():
		muonhitqdc=muonhits[i][plane]-shiftdict[int(muonhitsbar[i][plane])+int(plane)*10] if (args.shifted and (muonhitsbar[i][plane])+int(plane)*10 in shiftdict.keys()) else muonhits[i][plane]
		if args.shifted and (not (muonhitsbar[i][plane])+int(plane)*10 in shiftdict.keys()):
			continue
		p1func.Fill(muonhits[i][plane]-shiftdict[(muonhitsbar[i][plane])+int(plane)*10])
for i in range(int(len(muonhits)/2)):
	#only loop until half
	#pretend 2 of the events happened at the same time
	mockqdc={}
	for plane in muonhits[2*i].keys():	
		if plane in muonhits[2*i+1].keys():	
			dimuonhitqdc=muonhits[2*i][plane]+muonhits[2*i+1][plane]-shiftdict[(muonhitsbar[2*i][plane])+int(plane)*10]-shiftdict[(muonhitsbar[2*i+1][plane])+int(plane)*10] if (args.shifted and (muonhitsbar[2*i][plane])+int(plane)*10 in shiftdict.keys() and (muonhitsbar[2*i+1][plane])+int(plane)*10 in shiftdict.keys()) else muonhits[2*i][plane]+muonhits[2*i+1][plane]
			if args.shifted and (not ((muonhitsbar[2*i][plane])+int(plane)*10 in shiftdict.keys()) and ((muonhitsbar[2*i+1][plane])+int(plane)*10 in shiftdict.keys())):
				continue
			p2func.Fill(dimuonhitqdc)
print(p1func.GetEntries())
p1func.Scale(1./p1func.Integral())
p2func.Scale(1./p2func.Integral())
with ROOT.TFile('{}.root'.format(args.outfile),'RECREATE') as ofile:
	ofile.WriteObject(p1func,'p1')
	ofile.WriteObject(p2func,'p2')
print('process ended successfully')
