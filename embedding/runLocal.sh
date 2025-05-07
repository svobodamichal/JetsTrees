#!/bin/bash

export FASTJET='/gpfs/mnt/gpfs01/star/pwg/prozorov/install/fastjet-install'
export FASTJET_CONTRIB='/gpfs/mnt/gpfs01/star/pwg/prozorov/install/fjcontrib-1.046'
cons

# set pthat=(3 5 7 9 11 15 20 25 30 40 50 -1)
# #PYTHIA6 cross sections [mb] for pT-hat bins (3 5 7 9 11 15 20 25 30 40 50 -1) to be used as weights
# set weights=(1.616e+0 1.355e-01 2.288e-02 5.524e-03 2.203e-03 3.437e-04 4.681e-05 8.532e-06 2.178e-06 1.198e-07 6.939e-09)

pthatmin=50
pthatmax=-1
xweight=6.939e-09

root4star -l -b -q 'StRoot/macros/runPicoHFJetMaker.C("testPythia6picoDsts_pt50_-1.list","output_test",0,"PicoDst")'