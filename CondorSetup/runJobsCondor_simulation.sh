#!/bin/bash
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export COIN_FULL_INDIRECT_RENDERING=1
echo $VO_CMS_SW_DIR
source $VO_CMS_SW_DIR/cmsset_default.sh
cd /home/work/risharma/work/CMSSW_10_3_1/src
eval `scramv1 runtime -sh`
root -q -b -l /home/work/risharma/work/git/DisplacedMultileptons/CondorSetup/runana.C\(\"$1\",\"$2\",\"$3\",\"$4\",\"$5\",\"$6\"\,\"$7\"\)
