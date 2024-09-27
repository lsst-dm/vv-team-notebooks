#!/usr/bin/env bash

# Capture the start time.
start_time=$(date +%s)

welcome="OpsRehearsal4 Visit Summary Plot Generation"
underline=$(printf "%-${#welcome}s" | tr ' ' '-')
echo -e "$welcome\n$underline"

# I'm not sure why I'm having to explicitly set this, but loading the stack fails otherwise...
source $HOME/.bashrc
# echo "SHELL = ${SHELL}"
# echo "PATH= ${PATH}"
export SHELL="bash"
# echo "SHELL = ${SHELL}"

cmdLineOptions=$1

# Load the LSST environment if not already loaded.
if [ -x "$(command -v eups)" ]; then
    echo "eups is indeed setup..."
    LSST_DISTRIB=$(eups list | grep lsst_distrib | grep setup)
fi
if [ -n "$LSST_DISTRIB" ]; then
    echo -e "$LSST_DISTRIB\n"
else
    echo "Setting up stack with current weekly"
    # source /sdf/group/rubin/sw/w_latest/loadLSST.bash
    source /sdf/group/rubin/sw/w_latest/loadLSST.sh
    setup lsst_sitcom
    setup lsst_distrib
    LSST_DISTRIB=$(eups list | grep lsst_distrib | grep setup)
    if [ -n "$LSST_DISTRIB" ]; then
        echo -e "$LSST_DISTRIB\n"
    else
        echo "The lsst_distrib package is not set up. Exiting."
        exit 1
    fi
fi

cmd="/sdf/home/l/laurenma/Python/makeVisitSummaryPlots.py ${cmdLineOptions}"
echo "Making Plots running:"
echo ${cmd}
eval $cmd
echo "Done!"

# Capture the end time and calculate the total runtime.
end_time=$(date +%s)
runtime=$((end_time - start_time))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$((runtime % 60))
printf "\nTotal runtime: %02d:%02d:%02d\n" $hours $minutes $seconds
