import pandas as pd
import numpy as np
import ROOT
import sys

# Read the CSV file into a pandas DataFrame

#argv[0]:script name
#argv[1]:csvfile
#argv[2]:outputfile

if len(sys.argv)!=3:raise Exception("Wrong number of arguments: needs 2, given {}".format(len(sys.argv)-1))
if not sys.argv[1].endswith(".csv"):raise Exception("argument 1 must end in .csv")
if not sys.argv[2].endswith(".root"):raise Exception("argument 2 must end in .root")
    
df = pd.read_csv(sys.argv[1])

# Create a ROOT file to store the TTree
output_file = ROOT.TFile(sys.argv[2], 'RECREATE')

# Create a TTree
tree = ROOT.TTree('gst', 'gst')

# Create branches in the TTree for each variable
Ev = np.zeros(1, dtype=float)
pxv = np.zeros(1, dtype=float)
pyv = np.zeros(1, dtype=float)
pzv = np.zeros(1, dtype=float)
neu = np.zeros(1, dtype=int)
cc = np.zeros(1, dtype=bool)
nuel = np.zeros(1, dtype=bool)
vtxx = np.zeros(1, dtype=float)
vtxy = np.zeros(1, dtype=float)
vtxz = np.zeros(1, dtype=float)
vtxt = np.zeros(1, dtype=float)
El = np.zeros(1, dtype=float)
pxl = np.zeros(1, dtype=float)
pyl = np.zeros(1, dtype=float)
pzl = np.zeros(1, dtype=float)
Ef = np.zeros(2, dtype=float)
pxf = np.zeros(2, dtype=float)
pyf = np.zeros(2, dtype=float)
pzf = np.zeros(2, dtype=float)
nf = np.zeros(1, dtype=int)
pdgf = np.zeros(2, dtype=np.int32)

tree.Branch('Ev', Ev, 'Ev/D')
tree.Branch('pxv', pxv, 'pxv/D')
tree.Branch('pyv', pyv, 'pyv/D')
tree.Branch('pzv', pzv, 'pzv/D')
tree.Branch('neu', neu, 'neu/I')
tree.Branch('cc', cc, 'cc/O')
tree.Branch('nuel', nuel, 'nuel/O')
tree.Branch('vtxx', vtxx, 'vtxx/D')
tree.Branch('vtxy', vtxy, 'vtxy/D')
tree.Branch('vtxz', vtxz, 'vtxz/D')
tree.Branch('vtxt', vtxt, 'vtxt/D')
tree.Branch('El', El, 'El/D')
tree.Branch('pxl', pxl, 'pxl/D')
tree.Branch('pyl', pyl, 'pyl/D')
tree.Branch('pzl', pzl, 'pzl/D')
tree.Branch('Ef', Ef, 'Ef[2]/D')
tree.Branch('pxf', pxf, 'pxf[2]/D')
tree.Branch('pyf', pyf, 'pyf[2]/D')
tree.Branch('pzf', pzf, 'pzf[2]/D')
tree.Branch('nf', nf, 'nf/I')
tree.Branch('pdgf', pdgf, 'pdgf[2]/I')
opangle=[]
xlist=[]
ylist=[]
# Fill the TTree with the data from the DataFrame
for index, row in df[df['particle_type']==32].iterrows():
     #print("index {}\n row {} \n".format(index,row))
    #print(df.iloc[index+1])
    if df.shape[0]<index+2:
        continue
    rowmu=df.iloc[index+1]
    rowmu2=df.iloc[index+2]
    if rowmu['particle_type']==-rowmu2['particle_type']:
        Ev[0] = np.sqrt(row['px']**2+row['py']**2+row['pz']**2+row['m']**2)
        pxv[0] = row['px']
        pyv[0] = row['py']
        pzv[0] = row['pz']
        neu[0] = row['particle_type']
        cc[0] = False
        nuel[0] = False
        vtxx[0] = row['vx']/1000
        vtxy[0] = row['vy']/1000
        vtxz[0] = row['vz']/1000+2.8 #2.8 is the beginning of the detector
        vtxt[0] = row['vt']/1000+2.8
        El[0] = np.sqrt(row['px']**2+row['py']**2+row['pz']**2+row['m']**2)
        pxl[0] = row['px']
        pyl[0] = row['py']
        pzl[0] = row['pz']
        Ef[0] = np.sqrt(rowmu['px']**2+rowmu['py']**2+rowmu['pz']**2+rowmu['m']**2)
        Ef[1] = np.sqrt(rowmu2['px']**2+rowmu2['py']**2+rowmu2['pz']**2+rowmu2['m']**2)
        pxf[0] = rowmu['px']
        pxf[1] = rowmu2['px']
        pyf[0] = rowmu['py']
        pyf[1] = rowmu2['py']
        pzf[0] = rowmu['pz']
        pzf[1] = rowmu2['pz']
        nf[0] = 2
        pdgf[0]=rowmu['particle_type']
        pdgf[1]=rowmu2['particle_type']
        tree.Fill()
        opangle.append(np.arccos((rowmu['px']*rowmu2['px']+rowmu['py']*rowmu2['py']+rowmu['pz']*rowmu2['pz'])/(np.sqrt(rowmu['px']**2+rowmu['py']**2+rowmu['pz']**2)*np.sqrt(rowmu2['px']**2+rowmu2['py']**2+rowmu2['pz']**2))))
        xlist.append(row['vx']/1000)
        ylist.append(row['vy']/1000)
print("Tree fill concluded, saving to output file")  
# Write the TTree to the output file
output_file.Write()
output_file.Close()
hist = ROOT.TH1F("Opening Angle","Opening Angle for Dark Higgs decay to dimuon;angle (rad); event count", 100,min(opangle),max(opangle))
hxystart=ROOT.TH2D("XY-Vertex","Dark Higgs Decay Vertex XY-Plane;x(cm);y(cm)",39*5,-8,-39-8,39*5,15.5,15.5+39)
hxystart.SetMarkerStyle(2)
hxystart.SetMarkerSize(2)
for value in opangle:
    hist.Fill(value)
for i in range(len(xlist)):
    hxystart.Fill(xlist [i],ylist[i])
c1=ROOT.TCanvas()
hist.Draw()
c1.SaveAs("OpAngleFromCsv.png")
hxystart.Draw("COLZ")
c1.SaveAs("xydecayvertex.png")
