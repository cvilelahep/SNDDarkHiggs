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
parser.add_argument("--l",dest="legend",default="placeholder")
parser.add_argument("--o",dest="outfile",default="outfile")
parser.add_argument("--f",dest="digfile",nargs="*")
parser.add_argument("--dir",dest="rootdir",default=".")
parser.add_argument("--usingdir",action='store_true',dest='usingdir')
parser.add_argument("--w",action='store_true',dest='weighted')
parser.add_argument("--SciFicut",action='store_true',dest='scificut')
parser.add_argument("--UScut",action='store_true',dest='uscut')
parser.add_argument("--Vetocut",action='store_true',dest='vetocut')
parser.add_argument("--XYcut",action='store_true',dest='xycut')
parser.add_argument("--Anglecut",action='store_true',dest='anglecut')
parser.add_argument("--shift",action='store_true',dest='shifted')
args=parser.parse_args()
#dict of pfunction shifts per bar
shiftdict={}
shiftdict[1]=6
shiftdict[11]=4
shiftdict[21]=10
shiftdict[31]=-10
shiftdict[42]=50
shiftdict[3]=2
shiftdict[13]=20
shiftdict[23]=14
shiftdict[33]=-4
shiftdict[43]=26
shiftdict[6]=-4
shiftdict[16]=2
shiftdict[26]=32
shiftdict[36]=-12
shiftdict[46]=40
shiftdict[4]=-4
shiftdict[14]=26
shiftdict[24]=32
shiftdict[34]=-6
shiftdict[44]=46
shiftdict[2]=22
shiftdict[12]=12
shiftdict[22]=26
shiftdict[32]=-2
shiftdict[5]=4
shiftdict[15]=30
shiftdict[25]=30
shiftdict[35]=6
shiftdict[45]=28
shiftdict[47]=10
shiftdict[41]=26
shiftdict[27]=20
shiftdict[37]=-8
shiftdict[7]=20
shiftdict[17]=6
'''shiftdict[2]=-18
shiftdict[12]=-22
shiftdict[22]=-16
shiftdict[32]=-32
shiftdict[43]=2
shiftdict[0]=-18
shiftdict[10]=-22
shiftdict[20]=-22
shiftdict[30]=-30
shiftdict[40]=-4
shiftdict[6]=-26
shiftdict[16]=-24
shiftdict[26]=-22
shiftdict[36]=-32
shiftdict[46]=22
shiftdict[8]=-22
shiftdict[18]=-26
shiftdict[28]=-20
shiftdict[38]=-36
shiftdict[48]=6
shiftdict[5]=-20
shiftdict[15]=-30
shiftdict[25]=-20
shiftdict[35]=-34
shiftdict[45]=4
shiftdict[42]=26
shiftdict[3]=-18
shiftdict[13]=-26
shiftdict[23]=-22
shiftdict[33]=-34
shiftdict[44]=24
shiftdict[9]=-18
shiftdict[19]=-26
shiftdict[29]=-22
shiftdict[39]=-34
shiftdict[49]=10
shiftdict[1]=-18
shiftdict[11]=-24
shiftdict[21]=-20
shiftdict[31]=-32
shiftdict[41]=4
shiftdict[4]=-24
shiftdict[14]=-26
shiftdict[24]=-20
shiftdict[34]=-34
shiftdict[7]=-20
shiftdict[17]=-22
shiftdict[27]=-16
shiftdict[37]=-30
shiftdict[47]=-14
'''

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
	if realdataaproved and (not args.vetocut or not veto or random.random()<2.5*10**-6):
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
p1func={}
p2func={}
for i in range(len(muonhits)):
	num=0
	for plane in muonhits[i].keys():
		muonhitqdc=muonhits[i][plane]-shiftdict[int(muonhitsbar[i][plane])+int(plane)*10] if (args.shifted and (muonhitsbar[i][plane])+int(plane)*10 in shiftdict.keys()) else muonhits[i][plane]
		if args.shifted and not((muonhitsbar[i][plane])+int(plane)*10 in shiftdict.keys()):
			continue
		if int(muonhitsbar[i][plane])+int(plane)*10 in p1func.keys():
			p1func[int(muonhitsbar[i][plane])+int(plane)*10 ].Fill(muonhitqdc)
		else:
			p1func[int(muonhitsbar[i][plane])+int(plane)*10 ]=ROOT.TH1D("p1func_plane{}_bar{}","p1func_plane{}_bar{};hit_qdc;\%".format(plane,muonhitsbar[i][plane],plane,muonhitsbar[i][plane]),125,0,250)
			p1func[int(muonhitsbar[i][plane])+int(plane)*10 ].Fill(muonhitqdc)
with ROOT.TFile('{}.root'.format(args.outfile),'RECREATE') as ofile:
	for key in p1func.keys():
		if p1func[key].GetEntries()>100:
			p1func[key].Scale(1./p1func[key].Integral())
			ofile.WriteObject(p1func[key],'p1_plane{}_bar{}'.format(key/10,key%10))
print('process ended successfully')
