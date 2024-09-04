python /afs/cern.ch/work/h/hsantos/private/SNDThesis/DH/notrackselection.py --f $@ --o outfile --realdata
python $SNDSW_ROOT/shipLHC/run_muonRecoSND.py -f outfile.root -g /eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V4_2023.root -c passing_mu_DS -sc 1 -s ./ -hf linearSlopeIntercept -o
python /afs/cern.ch/work/h/hsantos/private/SNDThesis/DH/ABCDmethod.py --l "real_muon_events" --f outfile.root --realdata --o ABCDout --pfunc pfunc_allcuts.root --XYcut --Anglecut --UScut --SciFicut --shift
