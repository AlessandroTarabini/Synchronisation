#!/bin/bash

declare -a files=("\/store\/mc\/RunIISummer20UL17MiniAODv2\/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8\/MINIAODSIM\/106X_mc2017_realistic_v9-v2\/130000\/638CB4D1-A901-FE45-8AFE-21298AACAD4D.root"
"\/store\/mc\/RunIISummer20UL17MiniAODv2\/VBF_HToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8\/MINIAODSIM\/106X_mc2017_realistic_v9-v2\/2500000\/2B2874CE-8DC5-8A41-B4FB-B0B68AA02C7B.root"
"\/store\/mc\/RunIISummer20UL17MiniAODv2\/WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8\/MINIAODSIM\/106X_mc2017_realistic_v9-v2\/40000\/8FF12A6D-44D8-9145-90FC-12F8364F3324.root"
"\/store\/mc\/RunIISummer20UL17MiniAODv2\/WminusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8\/MINIAODSIM\/106X_mc2017_realistic_v9-v2\/250000\/EE928E68-C8F9-D045-BFC3-91AA32BEB8F4.root"
"\/store\/mc\/RunIISummer20UL17MiniAODv2\/ZH_HToZZ_4LFilter_M125_TuneCP5_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8\/MINIAODSIM\/106X_mc2017_realistic_v9-v2\/260000\/29EA6E9B-C0F2-E140-B3B0-C2E5E0FC36DA.root")

declare -a prodModes=("ggH_2017"
"VBF_2017"
"Wplus_2017"
"Wminus_2017"
"ZH_2017")

voms-proxy-init -voms cms
cp /tmp/x509up_u133688 /afs/cern.ch/user/a/atarabin/

for i in "${!files[@]}"; do
  if [ -d ${prodModes[$i]} ]; then rm -rf ${prodModes[$i]}; fi
  mkdir ${prodModes[$i]}
  cp analyzer_Sync_2017.py ${prodModes[$i]}
  cp batchScript.sh ${prodModes[$i]}
  cp condor.sub ${prodModes[$i]}
  cd ${prodModes[$i]}
  mkdir output error log
  sed -i "s/FOLDER/${prodModes[$i]}/g" condor.sub
  sed -i "s/MINIAOD_NAME/${files[$i]}/g" analyzer_Sync_2017.py
  condor_submit condor.sub
  cd ..
done
