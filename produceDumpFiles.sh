#!/bin/bash

declare -a prodModes=(
# "ggH_2016"
# "VBF_2016"
# "Wplus_2016"
# "Wminus_2016"
# "ZH_2016"
"ggH_2018"
"VBF_2018"
"Wplus_2018"
"Wminus_2018"
"ZH_2018"
"ttH_2018")

for i in "${!prodModes[@]}"; do
  python prepareSync.py -i ${prodModes[$i]}/ZZ4lAnalysis.root
  year=($(echo ${prodModes[$i]} | rev | cut -d"_" -f1  | rev))
  pMode=($(echo ${prodModes[$i]} | rev | cut -d"_" -f2-  | rev))
  pMode="$pMode.txt"
  cp CJLSTevents.txt /eos/user/a/atarabin/www/fiducial/synchUL/$year/$pMode
  cp CJLSTevents.txt ${prodModes[$i]}/.
done
