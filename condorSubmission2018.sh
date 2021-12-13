#!/bin/bash

declare -a files=(
"\/store\/mc\/RunIISummer20UL18MiniAODv2\/WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8\/MINIAODSIM\/106X_upgrade2018_realistic_v16_L1v1-v1\/270000\/01072974-E38F-FD48-A5DC-B30DE99991E0.root"
"\/store\/mc\/RunIISummer20UL18MiniAODv2\/ZH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8\/MINIAODSIM\/106X_upgrade2018_realistic_v16_L1v1-v2\/260000\/36985B2E-3FA2-7E4B-A5FC-D20D8B713766.root"
"\/store\/mc\/RunIISummer20UL18MiniAODv2\/ttH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8\/MINIAODSIM\/106X_upgrade2018_realistic_v16_L1v1-v2\/260000\/C0DA604A-0549-7844-B58B-846A0F0D9EDD.root"
)

declare -a prodModes=(
"Wplus_2018"
"ZH_2018"
"ttH_2018"
)

voms-proxy-init -voms cms
cp /tmp/x509up_u133688 /afs/cern.ch/user/a/atarabin/

for i in "${!files[@]}"; do
  if [ -d ${prodModes[$i]} ]; then rm -rf ${prodModes[$i]}; fi
  mkdir ${prodModes[$i]}
  cp analyzer_Sync_2018.py ${prodModes[$i]}
  cp batchScript.sh ${prodModes[$i]}
  cp condor.sub ${prodModes[$i]}
  cd ${prodModes[$i]}
  mkdir output log error
  sed -i "s/FOLDER/${prodModes[$i]}/g" condor.sub
  sed -i "s/2017/2018/g" batchScript.sh
  sed -i "s/MINIAOD_NAME/${files[$i]}/g" analyzer_Sync_2018.py
  condor_submit condor.sub
  cd ..
done
