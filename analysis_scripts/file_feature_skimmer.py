import ROOT
import numpy as np
import sys
import os
import random
import SndlhcGeo
from tqdm import tqdm
from statistics import mode
from argparse import ArgumentParser
from collections import defaultdict
ROOT.gStyle.SetOptStat(0)
parser=ArgumentParser()
parser.add_argument("--realdata",action='store_true',dest='realdata')
parser.add_argument("--l",dest="legend",default="placeholder")
parser.add_argument("--f",dest="digfile",nargs="*")
parser.add_argument("--o",dest="outfile",default="default")
parser.add_argument("--SciFicut",action='store_true',dest='scificut')
parser.add_argument("--UScut",action='store_true',dest='uscut')
parser.add_argument("--Anglecut",action='store_true',dest='anglecut')
parser.add_argument("--XYcut",action='store_true',dest='xycut')
parser.add_argument("--dir",dest="rootdir",default=".")
parser.add_argument("--usingdir",action='store_true',dest='usingdir')
parser.add_argument("--shift",action='store_true',dest='shifted')
parser.add_argument("--pfunc",dest="pfile",required=True)
args=parser.parse_args()
if not args.outfile:
    print('no outfile')
    exit()
    
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
	qdcthreshold=100
else:
	treeinfoname="cbmsim"
	qdcthreshold=80
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
#creating TTree
output_file = ROOT.TFile((f'{args.outfile}.root' if args.outfile!='default' else args.digfile[0].replace('.root','_skimmed.root')), 'RECREATE')
otree = ROOT.TTree('selection_features', 'selection_features')
#EVENT ID
branch_eventid= np.zeros(1, dtype=np.int32)
#SCIFI CUT
branch_maxdisttotrack= np.zeros(1, dtype=float)
branch_SciFihitsV= np.zeros(5, dtype=np.int32)
branch_SciFihitsH= np.zeros(5, dtype=np.int32)
#US CUT
branch_UShits= np.zeros(5, dtype=np.int32)
#XY and Angle cut
branch_veto=np.zeros(1, dtype=bool)
branch_vetocoords= np.zeros(2, dtype=float)
branch_polarangle= np.zeros(1, dtype=float)
#NLLR
branch_NLLR=np.zeros(1, dtype=float)

#Create Branches
otree.Branch('Event ID', branch_eventid, 'ID/I')
otree.Branch('Veto', branch_veto, 'veto/O')
otree.Branch('Max_Dist_to_Track', branch_maxdisttotrack, 'maxdist/D')
otree.Branch('V_SciFi_Hits', branch_SciFihitsV, 'vscifihits[5]/I')
otree.Branch('H_SciFi_Hits', branch_SciFihitsH, 'hscifihits[5]/I')
otree.Branch('US_Hits', branch_UShits, 'ushits[5]/I')
otree.Branch('Particle_XY_Veto', branch_vetocoords, 'vetoxy[2]/D')
otree.Branch('Polar_Angle', branch_polarangle, 'polarangle/D')
otree.Branch('NLLR', branch_NLLR, 'NLLR/D')
#now to fill the tree
for i_event, [eventdigi,eventreco] in tqdm(enumerate(zip(treedigi,treereco))) :
	realdataaproved=True
	a=ROOT.TVector3()
	b=ROOT.TVector3()
	veto=False
	if args.realdata:
		if eventdigi.EventHeader.GetBeamMode()!=11 or eventdigi.EventHeader.isIP1()!=True or eventdigi.EventHeader.GetEventTime()-prevtime<100:
			realdataaproved=False
	prevtime=eventdigi.EventHeader.GetEventTime()
	branch_veto[0]=False
	for hit in eventdigi.Digi_MuFilterHits:
		if hit.GetSystem()==1:
			branch_veto[0]=True
			veto=True
			break

	if realdataaproved:# and (not veto or random.random()<2.5*10**-6):
		if eventreco.Reco_MuonTracks.GetEntries()>0:
			for track in eventreco.Reco_MuonTracks:
				reject=False
				hitcountscifi=0
				hitcountx=0
				hitcounty=0
				passflag=True
				maxhittrackdiffx=0
				maxhittrackdiffy=0
				tempmufiupplane={1:0,2:0,3:0,4:0,0:0}
				tempmufidwplane={1:0,2:0,3:0,0:0}
				tempqdcmindw=[]
				tempqdcminup=[]
				associatedbar={}
				bestqdcperplaneup={} 
				minhittrackdiffperplaneup={}
				bestqdcperplanedw={} 
				minhittrackdiffperplanedw={}
				IDofminhitup={}
				IDofminhitdw={}
				distofbestqdcup={}
						
				for hit in eventdigi.Digi_MuFilterHits:
					if not hit.isValid():
						continue
					if hit.GetSystem()==3:
						tempmufidwplane[hit.GetPlane()]+=1
					if hit.GetSystem()==2:
						tempmufiupplane[hit.GetPlane()]+=1
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
				if args.uscut:
					if tempmufiupplane:
						if sum(tempmufiupplane.values())>5 or max(tempmufiupplane.values())>2:
							reject=True
							continue
					elif not tempmufiupplane:
						reject=True
						continue					
				station_stats = {}  # To store info about visited stations	
				for i in {1,2,3,4,5}:
						station_stats[i] = {
							"sum_x": 0,
							"sum_y": 0,
							"sum_z": 0,	
							"countx": 0,
							"county": 0,
							"countz": 0,
							"has_vertical": False,
							"has_horizontal": False,
						}
				scifipassx,scifipassy=False,False
				for hit in eventdigi.Digi_ScifiHits:						
					if not hit.isValid():
						continue							
					scifi.GetSiPMPosition(hit.GetChannelID(),a,b)
					# Track station stats
					station=hit.GetStation()
					if hit.isVertical():
						if station==5:
							scifipassx=True
						maxhittrackdiffx=max(maxhittrackdiffx,(getTrackAtZ(track,a.Z()).X()-(a.X()/2+b.X()/2))**2)						
						station_stats[station]["sum_x"]+=a.X()/2+b.X()/2
						station_stats[station]["countx"] += 1
						station_stats[station]["sum_z"]+=a.Z()/2+b.Z()/2
						station_stats[station]["countz"] += 1
						station_stats[station]["has_vertical"] = True
						
					else:
						if station==5:
							scifipassy=True
						maxhittrackdiffy=max(maxhittrackdiffy,(getTrackAtZ(track,a.Z()).Y()-(a.Y()/2+b.Y()/2))**2)			
						station_stats[station]["sum_y"]+=a.Y()/2+b.Y()/2
						station_stats[station]["county"] += 1
						station_stats[station]["sum_z"]+=a.Z()/2+b.Z()/2
						station_stats[station]["countz"] += 1
						station_stats[station]["has_horizontal"] = True	
				if args.scificut and not(scifipassy and scifipassx and max(maxhittrackdiffx,maxhittrackdiffy)<5):
					reject=True
					continue
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
				if args.xycut and not (np.logical_and(-44.5<vetocoords.X()<-10.5,18<vetocoords.Y()<52)):  
					reject=True
					continue
				if args.anglecut and not ((trackangleXZ*1000)**2+(trackangleYZ*1000)**2<=10):
					reject=True
					continue
				if not reject and station_stats and tempmufiupplane:
					branch_maxdisttotrack[0]=(np.sqrt(max(maxhittrackdiffx,maxhittrackdiffy)))
					branch_eventid[0]=eventdigi.EventHeader.GetEventNumber()
					branch_vetocoords[0]=(vetocoords.X())
					branch_vetocoords[1]=(vetocoords.Y())						
					if station_stats:
						branch_SciFihitsV[:]=[station_stats[i]["county"] for i in station_stats]
						branch_SciFihitsH[:]=[station_stats[i]["countx"] for i in station_stats]
					if tempmufiupplane:
						branch_UShits[:]=[tempmufiupplane[i] for i in tempmufiupplane]
					branch_polarangle[0]=np.sqrt(abs(trackangleXZ)**2+abs(trackangleYZ)**2)
					logcount=0
					for plane in bestqdcperplaneup.keys():
						muonhit_qdc=bestqdcperplaneup[plane]-shiftdict[int(associatedbar[plane])+int(plane)*10] if (args.shifted and (associatedbar[plane])+int(plane)*10 in shiftdict.keys()) else bestqdcperplaneup[plane]
						if args.shifted and not ((associatedbar[plane])+int(plane)*10 in shiftdict.keys()):
							continue
						logcount+=np.log(p1.GetBinContent(int(int(muonhit_qdc)/2)+1)) if p1.GetBinContent(int(int(muonhit_qdc)/2)+1)!=0 else 0
						logcount-=np.log(p2.GetBinContent(int(int(muonhit_qdc)/2)+1)) if p2.GetBinContent(int(int(muonhit_qdc)/2)+1)!=0 else 0
					branch_NLLR[0]=logcount
					otree.Fill()
output_file.Write()
output_file.Close()	

