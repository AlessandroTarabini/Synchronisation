#!/bin/bash
set -euo pipefail

if [ -z ${_CONDOR_SCRATCH_DIR+x} ]; then
  #running locally
  runninglocally=true
  _CONDOR_SCRATCH_DIR=$(mktemp -d)
  SUBMIT_DIR=$(pwd)
else
  runninglocally=false
  SUBMIT_DIR=$1
fi

cd $SUBMIT_DIR
eval $(scram ru -sh)

cp analyzer_Sync_2017.py $_CONDOR_SCRATCH_DIR
cd $_CONDOR_SCRATCH_DIR

pwd

echo 'Running at:' $(date)

cmsRunStatus=   #default for successful completion is an empty file
cmsRun analyzer_Sync_2017.py |& grep -v -e 'MINUIT WARNING' -e 'Second derivative zero' -e 'Negative diagonal element' -e 'added to diagonal of error matrix' > log.txt || cmsRunStatus=$?

echo -n $cmsRunStatus > exitStatus.txt
echo 'cmsRun done at: ' $(date) with exit status: ${cmsRunStatus+0}
gzip log.txt

export ROOT_HIST=0
if [ -s ZZ4lAnalysis.root ]; then
 root -q -b '${CMSSW_BASE}/src/ZZAnalysis/AnalysisStep/test/prod/rootFileIntegrity.r("ZZ4lAnalysis.root")'
else
 echo moving empty file
 mv ZZ4lAnalysis.root ZZ4lAnalysis.root.empty
fi

#delete mela stuff and $USER.cc
#I have no idea what $USER.cc is
rm -f br.sm1 br.sm2 ffwarn.dat input.DAT process.DAT "$USER.cc"

echo '...done at' $(date)

#note cping back is handled automatically by condor
if $runninglocally; then
  cp ZZ4lAnalysis.root* *.txt *.gz $SUBMIT_DIR
fi

exit $cmsRunStatus
