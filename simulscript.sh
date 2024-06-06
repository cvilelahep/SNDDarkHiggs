#must be called with the form
# source script.sh example.root so $i is the root file name
#1: simulation 

python $SNDSW_ROOT/shipLHC/run_simSND.py --Genie 4 -f $1 -i 0 -n 1000
#might give error, rerun
#python /eos/user/h/hsantos/SND/sndsw/shipLHC/run_simSND.py --Genie 4 -f $1 -i 0 -n 1000

#2 digitize

python $SNDSW_ROOT/shipLHC/run_digiSND.py -g geofile_full.Genie-TGeant4.root -f sndLHC.Genie-TGeant4.root

echo "loopEvents(save=True, auto=True, start=0, withHoughTrack = 3, withTrack = -1)" | python -i $SNDSW_ROOT/shipLHC/scripts/2dEventDisplay.py -p $(pwd)/ -f sndLHC.Genie-TGeant4_dig.root -g geofile_full.Genie-TGeant4.root

python $SNDSW_ROOT/shipLHC/run_muonRecoSND.py -f <FILE NAME> -g /eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V4_2023.root -c passing_mu_DS -sc 1 -s ./ -hf linearSlopeIntercept -o

