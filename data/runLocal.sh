#!/bin/bash

export FASTJET='/gpfs/mnt/gpfs01/star/pwg/prozorov/install/fastjet-install'
export FASTJET_CONTRIB='/gpfs/mnt/gpfs01/star/pwg/prozorov/install/fjcontrib-1.046'
cons
root4star -l -b -q 'StRoot/macros/runPicoHFJetMaker.C("./filelists/testPythia6picoDsts_pt50_-1_1file.list","test")'