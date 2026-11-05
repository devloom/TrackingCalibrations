#!/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
export MYINSTALL=/sphenix/u/dloomis/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

#printenv

echo "Done with setup"

runnumber=$1
prod=$2
type=$3
outDir=$4
inDir=$5

# runnumber: 78242
# prod: new_nocdbtag_v001
# type: run3auau
# outDir: Feb24_2026_LaminationFits_beamOffFieldOff
# inDir: Feb24_2026_LaserClusters_beamOffFieldOff

runShort=$(($1/100))


if [ ! -d "/sphenix/tg/tg01/hf/dloomis/CMDistortions/$outDir" ]
then
    mkdir -p /sphenix/tg/tg01/hf/dloomis/CMDistortions/$outDir
fi

filelist=run${runnumber}.list

rm $filelist

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`

for f in `ls /sphenix/tg/tg01/hf/dloomis/CMDistortions/${inDir}/DST_LASER_CLUSTER_${type}_${prod}-000${runnumber}-?????.root`
do
  echo $f >> $filelist
  if [[ `wc -l < $filelist` -gt 100 ]]
  then
    break;
  fi
done



# if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]
# then
#   cd $_CONDOR_SCRATCH_DIR
#   rsync -av $this_dir/* .
# else
#   echo condor scratch NOT set
#   exit -1
#  fi


# root -b -q Fun4All_LaminationFitting.C\(-1,\"${filelist}\",\"/sphenix/tg/tg01/hf/dloomis/CMDistortions/$outDir/\",\"$prod\",\"$type\"\)
root -b -q Fun4All_LaminationFitting.C\(-1,\"${filelist}\",\"/sphenix/tg/tg01/hf/dloomis/CMDistortions/$outDir/\",\"$prod\",\"$type\"\)
# root -b -q Fun4All_LaminationFitting.C\(-1,\"${filelist}\",\"./\",\"$prod\",\"$type\"\)

# cp Laminations_${type}_${prod}-000${runnumber}.root /sphenix/tg/tg01/hf/dloomis/CMDistortions/$outDir/
# cd $this_dir

rm $filelist
