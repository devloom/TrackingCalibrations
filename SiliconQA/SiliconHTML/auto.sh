#!/bin/bash

original_dir=$(pwd)

while [ 1 ]; do 
    # Start the execution timer
    SECONDS=0

    # clear the screen
    clear
    start_time=$(date "+%Y-%m-%d %H:%M:%S %Z")
    echo "The current execution started at $start_time"

    # Run the desired task
    root -l -q -b runsummary.C

    # Copy the status and results
    cp runstatus.html /sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/
    # cp ClusPhiChecker/* /sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/ClusPhi/
    # cp ChipOccupancyChecker/* /sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/ChipOccupancy/

    # Process the date-stamped Graph files
    current_time=$(date "+%Y%m%d")
    rm /sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/*Graph*
    cp *Graph.p* /sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary
    cd /sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary
    for file in *Graph.p*; do
        mv "$file" "${file//Graph/Graph_${current_time}}"
    done

    # Create the PDF summary
    gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=SummaryPlots.pdf -dBATCH *Graph_${current_time}.pdf

    cd $original_dir

    clear

    # Calculate execution time
    exec_time=$SECONDS
    echo "The script took $exec_time seconds to execute."

    # Calculate the next execution time (24 hours = 86400 seconds)
    next_execution_time=$((86400 - exec_time))

    echo "Sleeping for $next_execution_time seconds; will resume at $(date -d "+$next_execution_time seconds" "+%Y-%m-%d %H:%M:%S %Z")"

    # Sleep until the next 24-hour cycle
    sleep $next_execution_time
done
