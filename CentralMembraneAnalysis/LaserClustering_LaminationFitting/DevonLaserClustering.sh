#!/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
export MYINSTALL=/sphenix/user/dloomis/install/

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

#printenv

echo "Done with setup"

runnumber=$1
segment=$2
prod=$3
type=$4
dir=$5


runShort=$(($1/100))


if [ ! -d "/sphenix/tg/tg01/hf/dloomis/CMDistortions/$dir" ]
then
    mkdir -p /sphenix/tg/tg01/hf/dloomis/CMDistortions/$dir
fi

filelist=run${runnumber}_seg${segment}.list

for i in {00..23}; do
    # ls /sphenix/lustre01/sphnxpro/production/run3pp/physics/new_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_0/run_00079500_00079600/DST_STREAMING_EVENT_ebdc${i}_0_run3pp_new_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist
    # ls /sphenix/lustre01/sphnxpro/production/run3pp/physics/new_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_1/run_00079500_00079600/DST_STREAMING_EVENT_ebdc${i}_1_run3pp_new_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist

    # ls /sphenix/lustre01/sphnxpro/production/run3pp/physics/ana532_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_0/run_00079500_00079600/DST_STREAMING_EVENT_ebdc${i}_0_run3pp_ana532_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist
    # ls /sphenix/lustre01/sphnxpro/production/run3pp/physics/ana532_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_1/run_00079500_00079600/DST_STREAMING_EVENT_ebdc${i}_1_run3pp_ana532_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist


        
    # ls /sphenix/lustre01/sphnxpro/production/run3pp/physics/ana532_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_0/run_00079500_00079600/DST_STREAMING_EVENT_ebdc${i}_0_run3pp_ana532_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist
    # ls /sphenix/lustre01/sphnxpro/production/run3pp/physics/ana532_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_1/run_00079500_00079600/DST_STREAMING_EVENT_ebdc${i}_1_run3pp_ana532_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist

    # ls /sphenix/lustre01/sphnxpro/production/run3pp/physics/ana532_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_0/run_00081000_00081100/DST_STREAMING_EVENT_ebdc${i}_0_run3pp_ana532_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist
    # ls /sphenix/lustre01/sphnxpro/production/run3pp/physics/ana532_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_1/run_00081000_00081100/DST_STREAMING_EVENT_ebdc${i}_1_run3pp_ana532_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist

    # ls /sphenix/lustre01/sphnxpro/production/run3auau/calib/new_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_0/run_00078200_00078300/DST_STREAMING_EVENT_ebdc${i}_0_run3calib_new_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist
    # ls /sphenix/lustre01/sphnxpro/production/run3auau/calib/new_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_1/run_00078200_00078300/DST_STREAMING_EVENT_ebdc${i}_1_run3calib_new_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist

    # ls /sphenix/lustre01/sphnxpro/production/run3auau/calib/new_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_0/run_00075100_00075200/DST_STREAMING_EVENT_ebdc${i}_0_run3calib_new_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist
    # ls /sphenix/lustre01/sphnxpro/production/run3auau/calib/new_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_1/run_00075100_00075200/DST_STREAMING_EVENT_ebdc${i}_1_run3calib_new_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist

    # ls /sphenix/lustre01/sphnxpro/production/run3pp/calib/new_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_0/run_00084200_00084300/DST_STREAMING_EVENT_ebdc${i}_0_run3calib_new_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist
    # ls /sphenix/lustre01/sphnxpro/production/run3pp/calib/new_nocdbtag_v001/DST_STREAMING_EVENT_ebdc${i}_1/run_00084200_00084300/DST_STREAMING_EVENT_ebdc${i}_1_run3calib_new_nocdbtag_v001-000${runnumber}-${segment}.root >> $filelist

    # ls /sphenix/lustre01/sphnxpro/production/run3calib/calib/ana556_2025p013_v001/DST_STREAMING_EVENT_ebdc${i}_0/run_00075000_00075100/DST_STREAMING_EVENT_ebdc${i}_0_run3calib_ana556_2025p013_v001-000${runnumber}-${segment}.root >> $filelist
    # ls /sphenix/lustre01/sphnxpro/production/run3calib/calib/ana556_2025p013_v001/DST_STREAMING_EVENT_ebdc${i}_1/run_00075000_00075100/DST_STREAMING_EVENT_ebdc${i}_1_run3calib_ana556_2025p013_v001-000${runnumber}-${segment}.root >> $filelist

    ls /sphenix/lustre01/sphnxpro/production/run3calib/calib/ana556_2025p013_v001/DST_STREAMING_EVENT_ebdc${i}_0/run_00083300_00083400/DST_STREAMING_EVENT_ebdc${i}_0_run3calib_ana556_2025p013_v001-000${runnumber}-${segment}.root >> $filelist
    ls /sphenix/lustre01/sphnxpro/production/run3calib/calib/ana556_2025p013_v001/DST_STREAMING_EVENT_ebdc${i}_1/run_00083300_00083400/DST_STREAMING_EVENT_ebdc${i}_1_run3calib_ana556_2025p013_v001-000${runnumber}-${segment}.root >> $filelist
done

#beam data
root -b -q Fun4All_LaminationClustering.C\(-1,\"${filelist}\",\"/sphenix/tg/tg01/hf/dloomis/CMDistortions/$dir/\",\"DST_LASER_CLUSTER\",\"$prod\",\"$type\"\)

rm $filelist
