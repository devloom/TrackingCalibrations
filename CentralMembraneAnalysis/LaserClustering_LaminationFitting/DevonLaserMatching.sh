#!/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}/macros/detectors/sPHENIX/
export MYINSTALL=/sphenix/user/bkimelman/sPHENIX/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

#printenv

echo "Done with setup"

runnumber=$1
prod=$2
type=$3
skipOutliers=$4
manualInterp=$5
dir=$6
distDir=$7
clusDir=$8

#echo $segment
runShort=$(($1/100))

#current just commented out for a test
if [ ! -d "/sphenix/tg/tg01/jets/bkimelman/BenProduction/$dir" ]; then
    mkdir -p /sphenix/tg/tg01/jets/bkimelman/BenProduction/$dir
fi

filelist=""

if [ ${skipOutliers} == true ]; then
    if [ ${manualInterp} == true ]; then
	filelist=run${runnumber}_laserMatching_skipOutliers_manualInterp.list
    else
	filelist=run${runnumber}_laserMatching_skipOutliers.list
    fi
else
    if [ ${manualInterp} == true ]; then
	filelist=run${runnumber}_laserMatching_manualInterp.list
    else
	filelist=run${runnumber}_laserMatching.list
    fi
fi

distortionFile=""

if [ -f ${filelist} ]
then
    rm ${filelist}
fi

#ls /sphenix/tg/tg01/jets/bkimelman/BenProduction/Aug5_2025_LaserClusters_DC/DST_LASER_CLUSTER_${type}_${prod}-$(printf "%08d" $runnumber)-?????.root >> $filelist
#ls /sphenix/tg/tg01/jets/bkimelman/BenProduction/Aug5_2025_LaserClusters_DC/DST_LASER_CLUSTER_${type}_${prod}-$(printf "%08d" $runnumber)-00000.root >> $filelist
#ls /sphenix/tg/tg01/jets/bkimelman/BenProduction/Aug7_2025_LaserClusters_DC/DST_LASER_CLUSTER_${type}_${prod}-$(printf "%08d" $runnumber)-?????.root >> $filelist

if [[ -n "$distDir" ]]
then
    if [[ $distDir == "CDB" ]]
    then
	distortionFile="CDB"
    else
	distortionFile=/sphenix/tg/tg01/jets/bkimelman/BenProduction/${distDir}/Laminations_${type}_${prod}-$(printf "%08d" $runnumber).root
    fi
fi

if [[ -n "$clusDir" ]]
then
    ls /sphenix/tg/tg01/jets/bkimelman/BenProduction/$clusDir/DST_LASER_CLUSTER_${type}_${prod}-$(printf "%08d" $runnumber)-?????.root >> $filelist
else
    if [[ $type =~ "run2" ]]
    then
	for f in `ls  /sphenix/lustre01/sphnxpro/production/${type}/physics/${prod}/DST_TRKR_CLUSTER/run_000${runShort}00_000$((runShort+1))00/dst/DST_TRKR_CLUSTER_${type}_${prod}-000${runnumber}-?????.root`; do
            echo $f  >> $filelist
            if [[ `wc -l < $filelist` -gt 20 ]]
            then
		break;
            fi
	done
    elif [[ $type == "run3calib" ]]
    then
	for f in `ls  /sphenix/lustre01/sphnxpro/production/run3auau/calib/${prod}/DST_TRKR_CLUSTER/run_000${runShort}00_000$((runShort+1))00/DST_TRKR_CLUSTER_${type}_${prod}-000${runnumber}-?????.root`; do
            echo $f  >> $filelist
            if [[ `wc -l < $filelist` -gt 20 ]]
            then
		break;
            fi
	done
    elif [[ $type == "run3auau" || $type == "run3pp" ]]
    then
	for f in `ls  /sphenix/lustre01/sphnxpro/production/${type}/physics/${prod}/DST_TRKR_CLUSTER/run_000${runShort}00_000$((runShort+1))00/DST_TRKR_CLUSTER_${type}_${prod}-000${runnumber}-?????.root`; do
            echo $f  >> $filelist
            if [[ `wc -l < $filelist` -gt 200 ]]
            then
		break;
            fi
	done
    fi
fi

cat $filelist

root -b -q Fun4All_LaserMatching.C\(-1,\"${filelist}\",${skipOutliers},${manualInterp},\"${distortionFile}\",\"/sphenix/tg/tg01/jets/bkimelman/BenProduction/$dir/\"\)

rm $filelist
