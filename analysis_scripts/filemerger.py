import ROOT
import sys
import os
import SndlhcGeo
from tqdm import tqdm
from argparse import ArgumentParser
parser=ArgumentParser()
parser.add_argument("--realdata",action='store_true',dest='realdata')
parser.add_argument("--o",dest="outfile",default="outfile")
parser.add_argument("--f",dest="digfile",nargs="*")
parser.add_argument("--dir",dest="rootdir",default=".")
parser.add_argument("--usingdir",action='store_true',dest='usingdir')

args=parser.parse_args()


addlist=[]#This addlist will keep the event numbers of events to keep

#import files
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

#Get geofile in current dir
geoname=[filename for filename in os.listdir('.') if filename.startswith("geofile")]
geo=SndlhcGeo.GeoInterface(geoname[0])	
scifi=geo.modules['Scifi']		

prevtime=0#used for run data event pre-processing
for i_event, eventdigi in tqdm(enumerate(treedigi)) :
	realdataaproved=True
	#check whether the run data event passes pre-processing
	if args.realdata:
		if eventdigi.EventHeader.GetBeamMode()!=11 or eventdigi.EventHeader.isIP1()!=True or eventdigi.EventHeader.GetEventTime()-prevtime<100:
			realdataaproved=False
	prevtime=eventdigi.EventHeader.GetEventTime()
	if realdataaproved:
		addlist.append(i_event)
print("adding {} events to new file".format(len(addlist)))

#create new file by copying an empty version of the old file structure and then filling with approved events
filedigi=treedigi.GetFile()
with ROOT.TFile("{}.root".format(args.outfile),'RECREATE') as newfiledigi:
	for key in filedigi.GetListOfKeys():
		if key.GetName()!=treeinfoname and key.GetName()!="cbmsim":# to retain structure, copy all non-tree components of the file
			newfiledigi.WriteObject(filedigi.Get(key.GetName()),key.GetName())
	newtreedigi=treedigi.CloneTree(0)#clone empty version of tree
	#Fill with approved events
	for i in addlist:
		treedigi.GetEntry(i)
		newtreedigi.Fill()
	newfiledigi.WriteObject(newtreedigi,treeinfoname)
