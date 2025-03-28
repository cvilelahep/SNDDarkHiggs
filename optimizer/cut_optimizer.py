import ROOT
import numpy as np
from tqdm import tqdm
import tracemalloc
import time
import pickle
import scipy.stats
import scipy
import os
from skopt import gp_minimize
from skopt import dump, load
import random
import itertools
import gc
from argparse import ArgumentParser
import csv

#TREE STRUCTURE

#otree.Branch('Event ID', branch_eventid, 'ID/I')
#otree.Branch('Max_Dist_to_Track', branch_maxdisttotrack, 'maxdist/D')
#otree.Branch('V_SciFi_Hits', branch_SciFihitsV, 'vscifihits[5]/I')
#otree.Branch('H_SciFi_Hits', branch_SciFihitsH, 'hscifihits[5]/I')
#otree.Branch('US_Hits', branch_UShits, 'ushits[5]/I')
#otree.Branch('Particle_XY_Veto', branch_vetocoords, 'vetoxy[2]/D')
#otree.Branch('Polar_Angle', branch_polarangle, 'polarangle/D')
#otree.Branch('NLLR', branch_NLLR, 'NLLR/D')
MAX_FLOAT=65504

LEAF_NAME={'Event ID':'ID',
         'Max_Dist_to_Track': 'maxdist',
        'V_SciFi_Hits':'vscifihits',
        'H_SciFi_Hits':'hscifihits',
        'US_Hits':'ushits',
        'Particle_XY_Veto':'vetoxy',
        'Polar_Angle':'polarangle',
        'NLLR':'NLLR'}
POISSON_CDF_THRESHOLD={0:2.31,
                       1:3.89,
                       2:5.33,
                       3:6.69}

VETO_XLEFT=-47.0
VETO_XRIGHT=-8
VETO_YLEFT=15.5
VETO_YRIGHT=54.5
NU_FILESIZE=8100#Fb-1
MU_FILESIZE=0.79#Fb-1
ROOT.gStyle.SetOptStat(0)
parser=ArgumentParser()
parser.add_argument("--sigf",dest="sigfile",nargs="*")
parser.add_argument("--save", dest="saving", nargs="*", default=None)
parser.add_argument("--load",dest="loading", nargs="*", default=None)
parser.add_argument("--nubkgf",dest="nubkgfile",nargs="*")
parser.add_argument("--mubkgf",dest="mubkgfile",nargs="*")
parser.add_argument("--mubkgdir",dest="mubkgdir",default=None)
args=parser.parse_args()
#parsing save argument
if args.saving is None:
    args.saving = False  # Case: --save is NOT used at all
elif len(args.saving) == 0:
    args.saving = ["DHArrs", "NuArrs", "DataArrs"]  # Case: --save used without arguments
elif len(args.saving) != 3:
    parser.error("--save requires exactly 3 arguments.")
#parsing load argument    
if args.loading is None:
    args.loading = False  # Case: --load is NOT used at all
elif len(args.loading) == 0:
    args.loading = ["DHArrs.pkl", "NuArrs.pkl", "DataArrs.pkl"]  # Case: --load used without arguments
elif len(args.loading) != 3:
    parser.error("--load requires exactly 3 arguments.")
tracemalloc.start()
sigtree=ROOT.TChain("selection_features")
nubkgtree=ROOT.TChain("selection_features")
mubkgtree=ROOT.TChain("selection_features")
if args.sigfile:
    for file in args.sigfile:
        print(f"Adding {file} to DH signal chain ")
        sigtree.Add(file)
if args.nubkgfile:
    for file in args.nubkgfile:
        print(f"Adding {file} to neutrino background chain ")
        nubkgtree.Add(file)
if args.mubkgfile:
    for file in args.mubkgfile:
        print(f"Adding {file} to rundata background chain ")
        mubkgtree.Add(file)
if args.mubkgdir:
    for file in [filename for filename in os.listdir(args.mubkgdir[:-1]) if filename.endswith(".root")]:
    #for file in [filename for filename in os.listdir(args.mubkgdir[:-1]) if (match := re.search(r'(\d+)\.root$', filename)) and int(match.group(1)) % 4 == 0]:
        print(f"Adding {file} to rundata background chain ")
        mubkgtree.Add(f'{args.mubkgdir}{file}')
sigtree_stats={}
nubkgtree_stats={}
mubkgtree_stats={}
sorted_sigtree_stats={}
sorted_sigtree_index={}
sorted_nubkgtree_stats={}
sorted_nubkgtree_index={}
sorted_mubkgtree_stats={}
sorted_mubkgtree_index={}    

def append_event_to_csv(event_features, csv_file):
    # Check if the file exists to determine if headers need to be written
    file_exists = False
    try:
        with open(csv_file, 'r'):
            file_exists = True
    except FileNotFoundError:
        pass

    # Open the file in append mode
    with open(csv_file, mode='a', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=event_features.keys())

        # Write headers if the file is new
        if not file_exists:
            writer.writeheader()

        # Write the event features as a new row
        writer.writerow(event_features)
        
def get_threshold_metrics(tree,tree_stats:dict):
    #distance to track
    b_maxdist=[]
    #US hits
    b_maxus=[]
    b_totus=[]
    #SciFi hits
    b_scifi_lastplane=[]
    #Veto 
    b_vetoflag=[]
    b_vetox=[]
    b_vetoy=[]
    #Others
    b_polarangle=[]
    b_NLLR=[]    
    
    for i_event, event in tqdm(enumerate(tree)):
        #PRE-SELECTION
        if (
            event.maxdist>50 or
            max(event.ushits)>4 or                                                                                                                                                                                                                                                            
            sum(event.ushits)>20 or
            event.polarangle>100 #or
            #event.NLLR>5
            ):
            continue
        #distance to track
        b_maxdist.append(np.float16(np.sqrt(event.maxdist)))
        #with warnings.catch_warnings(record=True) as w:
        #    warnings.simplefilter("always", RuntimeWarning)  # Capture all RuntimeWarnings
        #    b_maxdist=np.float16(event.maxdist)
#
        #    # Check if an overflow warning occurred
        #    if any("overflow" in str(warning.message).lower() for warning in w):
        #        print(f"Overflow warning: Value too large for float16! value= {event.maxdist}")       
        #US hits
        b_maxus.append(max(event.ushits))
        b_totus.append(sum(event.ushits))
        #SciFi hits
        #b_maxscifi=max([event.vscifihits[i]+event.hscifihits[i] for i in range(5)])
        #b_totscifi=sum(event.vscifihits)+sum(event.hscifihits)
        b_scifi_lastplane.append(event.vscifihits[-1]+event.hscifihits[-1] if (event.vscifihits[-1]!=0 and event.hscifihits[-1]!=0) else 9999)
        #Veto     
        b_vetoflag.append(event.veto)
        b_vetox.append(event.vetoxy[0])
        b_vetoy.append(event.vetoxy[1])
        #Others
        b_polarangle.append(np.float16(np.sqrt(event.polarangle)))
        #with warnings.catch_warnings(record=True) as w:
        #    warnings.simplefilter("always", RuntimeWarning)  # Capture all RuntimeWarnings
        #    b_polarangle=np.float16(event.polarangle)
        #    # Check if an overflow warning occurred
        #    if any("overflow" in str(warning.message).lower() for warning in w):
        #        print(f"Overflow warning: Value too large for float16! value= {event.polarangle}")               
        b_NLLR.append(np.float16(event.NLLR))
    #distance to track
    tree_stats['b_maxdist']=np.array(b_maxdist,dtype=np.float16)
    #US hits
    tree_stats['b_maxus']=np.array(b_maxus, dtype=np.int8)
    tree_stats['b_totus']=np.array(b_totus, dtype=np.int8)
    #SciFi hits
    tree_stats['b_scifi_lastplane']=np.array(b_scifi_lastplane, dtype=np.int16)
    #Veto 
    tree_stats['b_vetoflag']=np.array(b_vetoflag,dtype=bool)
    tree_stats['b_vetox']=np.array(b_vetox,dtype=np.float16)
    tree_stats['b_vetoy']=np.array(b_vetoy,dtype=np.float16)
    #Others
    tree_stats['b_polarangle']=np.array(b_polarangle,dtype=np.float16)
    tree_stats['b_NLLR']=np.array(b_NLLR,dtype=np.float16)   
    print(f'Events selected: {len(tree_stats["b_maxdist"])}')
def sort_tree_stats(tree_stats,sorted_index,sorted_stats):    
    featlist=list(tree_stats.keys())
    for feature in featlist:
        sorted_index[feature]=np.argsort(tree_stats[feature])
        sorted_stats[feature]=tree_stats[feature][sorted_index[feature]]    
        del tree_stats[feature]
          
def var_plot(stat:str,n,min,max,legend,unit='',outfile='plots'):
    plots={}
    for name,tree_stats in [['signal',sorted_sigtree_stats],['nu',sorted_nubkgtree_stats],['mu',sorted_mubkgtree_stats]]:
        plots[name]=ROOT.TH1D(f"{legend}_{name}",f"{legend} for {name} events;{legend}{unit};count",n,min,max)
        for i in tree_stats[stat]:
            if i != 0:
                plots[name].Fill(i)
    with ROOT.TFile('{}.root'.format(outfile),'RECREATE') as ofile:
        for plot in plots.values():
            plot.Scale(1./plot.Integral())
            ofile.WriteObject(plot,legend)
def print_data_event(i,sorted_stats,sorted_index):
    event_features={}
    event_features['ID']=i
    for feature in sorted_stats.keys():
        event_pos = np.where(sorted_index[feature] == i)[0][0]
        event_features[feature]=sorted_stats[feature][event_pos]
    append_event_to_csv(event_features, 'event_features.csv')
def upperLimitFlatPrior(CL, n_obs, b):
    p = 1 - (1-CL)*(1 - scipy.stats.chi2.cdf(x = 2*b, df = 2*(n_obs+1)))
    return 1/2*scipy.stats.chi2.ppf(p, 2*(n_obs+1)) - b              
def check_thresholds(NLLRcut=-6,maxdist_cut=5,maxus_cut=2,totus_cut=5,scifi_lastplane_cut=20,polarangle_cut=10,x_offsetR=2.5,y_offsetR=2.5,x_offsetL=None,y_offsetL=None,printflag=False):
    if not x_offsetL:
        x_offsetL=x_offsetR
    if not y_offsetL:
        y_offsetL=y_offsetR 
    print(locals()) if printflag else None
    #Signal
    passing_sig={}
    #SciFi, US and Angle cut        
    for cut, threshold in [['b_maxdist',maxdist_cut],
                           ['b_scifi_lastplane',scifi_lastplane_cut],
                           ['b_maxus',maxus_cut],
                           ['b_totus',totus_cut],
                           ['b_polarangle',polarangle_cut],
                           ['b_NLLR',NLLRcut],                           
                           ['b_vetoflag',False]]:
        index=np.searchsorted(sorted_sigtree_stats[cut], threshold, side="right")
        passing_sig[cut] = set(sorted_sigtree_index[cut][:index])
    #XY cut
    #offset represents the frame cut additionally from the Veto dimensions
    idx_xright=np.searchsorted(sorted_sigtree_stats['b_vetox'], VETO_XRIGHT-x_offsetR, side="right")
    idx_xleft=np.searchsorted(sorted_sigtree_stats['b_vetox'], VETO_XLEFT+x_offsetL, side="left")
    idx_yright=np.searchsorted(sorted_sigtree_stats['b_vetoy'], VETO_YRIGHT-y_offsetR, side="right")
    idx_yleft=np.searchsorted(sorted_sigtree_stats['b_vetoy'], VETO_YLEFT+y_offsetL, side="left")
    passing_sig['b_vetox'] = set(sorted_sigtree_index['b_vetoy'][idx_xleft:idx_xright])
    passing_sig['b_vetoy'] = set(sorted_sigtree_index['b_vetoy'][idx_yleft:idx_yright])
    passing_all_sig=set.intersection(*passing_sig.values())
    sig_eff=len(passing_all_sig)/len(sorted_sigtree_stats['b_maxdist'])
    print(f"Signal efficiency: {sig_eff}") if printflag else None
    
    #Neutrino Background
    passing_nubkg={}
    passing_all_nubkg=None
    #SciFi, US and Angle cut        
    for cut, threshold in [['b_maxdist',maxdist_cut],
                           ['b_scifi_lastplane',scifi_lastplane_cut],
                           ['b_maxus',maxus_cut],
                           ['b_totus',totus_cut],
                           ['b_polarangle',polarangle_cut],
                           ['b_NLLR',NLLRcut],                           
                           ['b_vetoflag',False]
                           ]:
        passed_current_cut=set(sorted_mubkgtree_index[cut][:index])
        passing_all_nubkg=passed_current_cut if not passing_all_nubkg else set.intersection(passing_all_nubkg,passed_current_cut)
        del passed_current_cut
    #XY cut
    #offset represents the frame cut additionally from the Veto dimensions
    idx_xright=np.searchsorted(sorted_nubkgtree_stats['b_vetox'], VETO_XRIGHT-x_offsetR, side="right")
    idx_xleft=np.searchsorted(sorted_nubkgtree_stats['b_vetox'], VETO_XLEFT+x_offsetL, side="left")
    idx_yright=np.searchsorted(sorted_nubkgtree_stats['b_vetoy'], VETO_YRIGHT-y_offsetR, side="right")
    idx_yleft=np.searchsorted(sorted_nubkgtree_stats['b_vetoy'], VETO_YLEFT+y_offsetL, side="left")
    passing_nubkg['b_vetox'] = set(sorted_nubkgtree_index['b_vetoy'][idx_xleft:idx_xright])
    passing_nubkg['b_vetoy'] = set(sorted_nubkgtree_index['b_vetoy'][idx_yleft:idx_yright])
    passing_all_nubkg=set.intersection(passing_all_nubkg,passing_nubkg['b_vetox'],passing_nubkg['b_vetoy'])
    del passing_nubkg['b_vetox']
    del passing_nubkg['b_vetoy']
    nu_bkg=len(passing_all_nubkg)*200/NU_FILESIZE
    del passing_all_nubkg
    print(f'Nu background= {nu_bkg}') if printflag else None
    #print(f"Background efficiency: {len(passing_all_nubkg)/len(sorted_nubkgtree_stats['b_maxdist'])}, {len(passing_all_nubkg)} events survived")
    
    #Muon Background
    passing_mubkg={}
    passing_preABCD_mubkg=None  # Initializing variable
    
    for cut, threshold in [
        ['b_maxdist', maxdist_cut],
        ['b_scifi_lastplane', scifi_lastplane_cut],
        ['b_maxus', maxus_cut],
        ['b_totus', totus_cut],
        ['b_polarangle', polarangle_cut],
        ['b_NLLR', NLLRcut],
        ['b_vetoflag', False]
    ]:
        index = np.searchsorted(sorted_mubkgtree_stats[cut], threshold, side="right")
        if cut in ['b_NLLR', 'b_vetoflag']:
            passing_mubkg[cut] = set(sorted_mubkgtree_index[cut][:index])
        else:
            passed_current_cut = set(sorted_mubkgtree_index[cut][:index])
            if passing_preABCD_mubkg is None:
                passing_preABCD_mubkg = passed_current_cut
            else:
                passing_preABCD_mubkg.intersection_update(passed_current_cut)  # In-place intersection
            del passed_current_cut  # Free memory as soon as it's used

    # XY cut
    idx_xright = np.searchsorted(sorted_mubkgtree_stats['b_vetox'], VETO_XRIGHT - x_offsetR, side="right")
    idx_xleft = np.searchsorted(sorted_mubkgtree_stats['b_vetox'], VETO_XLEFT + x_offsetL, side="left")
    idx_yright = np.searchsorted(sorted_mubkgtree_stats['b_vetoy'], VETO_YRIGHT - y_offsetR, side="right")
    idx_yleft = np.searchsorted(sorted_mubkgtree_stats['b_vetoy'], VETO_YLEFT + y_offsetL, side="left")

    passing_mubkg['b_vetox'] = set(sorted_mubkgtree_index['b_vetox'][idx_xleft:idx_xright])
    passing_mubkg['b_vetoy'] = set(sorted_mubkgtree_index['b_vetoy'][idx_yleft:idx_yright])

    # Combine spatial cuts
    passing_preABCD_mubkg.intersection_update(passing_mubkg['b_vetox'])
    passing_preABCD_mubkg.intersection_update(passing_mubkg['b_vetoy'])

    # Clean up
    del passing_mubkg['b_vetox']
    del passing_mubkg['b_vetoy']
    #regionA=set([i for i in passing_preABCD_mubkg if (i in passing_mubkg['b_vetoflag'] and i in passing_mubkg['b_NLLR'])])
    regionB=len([i for i in passing_preABCD_mubkg if (i not in passing_mubkg['b_vetoflag'] and i in passing_mubkg['b_NLLR'])])
    regionC=len([i for i in passing_preABCD_mubkg if (i in passing_mubkg['b_vetoflag'] and i not in passing_mubkg['b_NLLR'])])
    regionD=len([i for i in passing_preABCD_mubkg if (i not in passing_mubkg['b_vetoflag'] and i not in passing_mubkg['b_NLLR'])])
    del passing_preABCD_mubkg
    del passing_mubkg
    if regionD==0:
        return 99999
    mu_bkg_expect=regionB*float(regionC)/regionD*200/MU_FILESIZE
    print(f'ABCD:B={regionB},C={regionC} and D={regionD}') if printflag else None
    print(f'Mu expected background={mu_bkg_expect}') if printflag else None
    if mu_bkg_expect>=5:
        return 99999
    optimization=0
    for i in range(15):
        optimization+=scipy.stats.poisson.pmf(i, mu_bkg_expect+nu_bkg)*upperLimitFlatPrior(CL=0.9,n_obs=i,b=mu_bkg_expect+nu_bkg)
    optimization=optimization/sig_eff if sig_eff!=0 else 99999
    print(f'{optimization} DH events needed for background exclusion at 90%') if printflag else None
    return optimization
 
         
    
    

    

if args.loading:
    print("Loading sorted data")
    start=time.time()
    gc.disable()
    with open(args.loading[0], "rb") as f:
        sorted_sigtree_index,sorted_sigtree_stats = pickle.load(f)
    print(f"Signal data loaded; {time.time()-start}")
    with open(args.loading[1], "rb") as f:
        sorted_nubkgtree_index,sorted_nubkgtree_stats = pickle.load(f)
    print(f"Neutrino data loaded; {time.time()-start}")
    with open(args.loading[2], "rb") as f:
        sorted_mubkgtree_index,sorted_mubkgtree_stats = pickle.load(f)
    print(f"Run data loaded; {time.time()-start}")
    end=time.time()
    gc.enable()
    print(f'Arrays loaded successfuly in {end-start} seconds')
else:
    print("getting signal data")
    get_threshold_metrics(sigtree,sigtree_stats)
    print("getting nu data")
    get_threshold_metrics(nubkgtree,nubkgtree_stats)
    print("getting mu data")
    get_threshold_metrics(mubkgtree,mubkgtree_stats)
    print("sorting signal data")
    sort_tree_stats(sigtree_stats,sorted_sigtree_index,sorted_sigtree_stats)
    print("sorting nu data")
    sort_tree_stats(nubkgtree_stats,sorted_nubkgtree_index,sorted_nubkgtree_stats)
    print("sorting mu data")
    sort_tree_stats(mubkgtree_stats,sorted_mubkgtree_index,sorted_mubkgtree_stats)
    if args.saving:
        print(f"Saving sorted stats: {args.saving[0]},{args.saving[1]},{args.saving[2]}")
        with open(f"{args.saving[0]}.pkl", "wb") as f:
            pickle.dump((sorted_sigtree_index,sorted_sigtree_stats), f)
        with open(f"{args.saving[1]}.pkl", "wb") as f:
            pickle.dump((sorted_nubkgtree_index,sorted_nubkgtree_stats), f)
        with open(f"{args.saving[2]}.pkl", "wb") as f:
            pickle.dump((sorted_mubkgtree_index,sorted_mubkgtree_stats), f)

         
            
GRIDSEARCH=False
if GRIDSEARCH:
    NLLR_cut_list=[-7,-6,-5]
    maxdist_cut_list=[np.sqrt(5),3,10]
    maxus_cut_list=[1,2,3]
    totus_cut_list=[5,6,7]
    scifi_lastplane_cut_list=[50,10]
    polarangle_cut_list=[5,10,100]
    x_offset_list=[2.5]
    y_offset_list=[2.5]
    best_combo=None
    best_opt=99999
    counter=0
    for threshold_combo in tqdm(itertools.product(NLLR_cut_list,maxdist_cut_list,maxus_cut_list,totus_cut_list,scifi_lastplane_cut_list,polarangle_cut_list,x_offset_list,y_offset_list),total=len([i for i in itertools.product(NLLR_cut_list,maxdist_cut_list,maxus_cut_list,totus_cut_list,scifi_lastplane_cut_list,polarangle_cut_list,x_offset_list,y_offset_list)]),desc="Processing..."):
        print("Beginning gridsearch") if counter==0 else None
        counter+=1
        optimization=check_thresholds(*threshold_combo)
        if optimization<best_opt:
            best_opt=optimization
            best_combo=threshold_combo
        if counter%50==0:
            print(counter,':',best_opt, '\n', best_combo)
    check_thresholds(*best_combo,printflag=True)

#Initialize with known good values
OPTIMIZER=True
if OPTIMIZER:
    OPTIMIZE_XY=False
    initial_input=[-7, 5, 1, 5, 60, 10, 1.5, 1.5] if OPTIMIZE_XY else [-7, 5, 1, 5, 60, 10]
    check_thresholds(*initial_input, x_offsetL=None, y_offsetL=None, printflag=True)
    bounds=[(-8.0,-4.0),(0.0,50.0),(1,5),(5,20),(0,100),(0.,100.),(0.,4),(0.,4)] if OPTIMIZE_XY else [(-8.0,-4.0),(0.0,50.0),(1,5),(5,20),(0,100),(0.,100.)]
    print("Beginning differential evolution of optimization metric")
    def objective_function(x):
        return check_thresholds(*x, x_offsetL=None, y_offsetL=None, printflag=False)    
    result = gp_minimize(objective_function,            # the function to minimize
                    dimensions=bounds,      # the bounds on each dimension of x
                    x0=initial_input,            # the starting point
                    acq_func="LCB",     # the acquisition function (optional)
                    n_calls=100,         # the number of evaluations of f including at x0
                    n_random_starts=4,  # the number of random initial points
                    random_state=2701,
                    verbose=True)          
    check_thresholds(*result.x, x_offsetL=None, y_offsetL=None, printflag=True)

# After the code runs, check memory usage
snapshot = tracemalloc.take_snapshot()

# Print the top 10 memory usage statistics
top_stats = snapshot.statistics('lineno')
for stat in top_stats[:10]:
    print(stat)




    

      




			
             
  






