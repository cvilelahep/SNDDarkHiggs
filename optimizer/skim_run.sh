python filemerger.py --f $@ --o outfile --realdata
python $SNDSW_ROOT/shipLHC/run_muonRecoSND.py -f outfile.root -g /eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V4_2023.root -c passing_mu_DS -sc 1 -s ./ -hf linearSlopeIntercept -o
python file_feature_skimmer.py --f outfile.root --pfunc pfunc_allcuts.root  --realdata --shift 