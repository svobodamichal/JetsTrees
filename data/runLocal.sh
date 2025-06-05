#!/bin/bash

# export FASTJET='/gpfs/mnt/gpfs01/star/pwg/prozorov/install/x64/fastjet-install'
# export FASTJET_CONTRIB='/gpfs/mnt/gpfs01/star/pwg/prozorov/install/x64/fjcontrib-1.046'

# export PATH=$FASTJET/bin:$PATH
# export LD_LIBRARY_PATH=$FASTJET/lib:$LD_LIBRARY_PATH
# export CPLUS_INCLUDE_PATH=$FASTJET/include:$CPLUS_INCLUDE_PATH

cons
root4star -l -b -q 'StRoot/macros/runPicoHFJetMaker.C("/star/data105/embedding/Run14_AuAu200_Pythia6_PicoDst/pt50_-1_0/out/st_physics_adc_15097061_raw_0000000.picoDst.root","test")'
