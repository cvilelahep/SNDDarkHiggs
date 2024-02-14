import pandas as pd
import numpy as np
import ROOT

# Read the CSV file into a pandas DataFrame
df = pd.read_csv('input.csv')

# Create a ROOT file to store the TTree
output_file = ROOT.TFile('output.root', 'RECREATE')

# Create a TTree
tree = ROOT.TTree('tree', 'Data')

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
pdgf = np.zeros(2, dtype=int)

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

# Fill the TTree with the data from the DataFrame
for index, row in df.iterrows():
    Ev[0] = row['Ev']
    pxv[0] = row['pxv']
    pyv[0] = row['pyv']
    pzv[0] = row['pzv']
    neu[0] = row['neu']
    cc[0] = row['cc']
    nuel[0] = row['nuel']
    vtxx[0] = row['vtxx']
    vtxy[0] = row['vtxy']
    vtxz[0] = row['vtxz']
    vtxt[0] = row['vtxt']
    El[0] = row['El']
    pxl[0] = row['pxl']
    pyl[0] = row['pyl']
    pzl[0] = row['pzl']
    Ef[0] = row['Ef'][0]
    Ef[1] = row['Ef'][1]
    pxf[0] = row['pxf'][0]
    pxf[1] = row['pxf'][1]
    pyf[0] = row['pyf'][0]
    pyf[1] = row['pyf'][1]
    pzf[0] = row['pzf'][0]
    pzf[1] = row['pzf'][1]
    nf[0] = row['nf']
    pdgf[0] = row['pdgf'][0]
    pdgf[1] = row['pdgf'][1]
    tree.Fill()

# Write the TTree to the output file
output_file.Write()
output_file.Close()
