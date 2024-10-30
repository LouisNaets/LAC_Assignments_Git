#!/bin/sh
#
#BSUB -J postprocess_hawc2
#BSUB -q hpc 
#BSUB -n 1
#BSUB -W 03:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2GB]"
#BSUB -o lsf_log/postprocess_hawc2_output_%J.out
#BSUB -e lsf_log/postprocess_hawc2_error_%J.err 

# CHECK BEFORE SUBMITTING THIS SCRIPT TO THE CLUSTER:
#   1. The filenames, directories, and settings in postprocess_hawc2.py are correct.
#   2. The results folders, the postprocess_hawc.py file, and this .sh are all in the same directory.
#   3. You have installed the lacbox as instructed during lecture.
#
# This script should save the results files as specified in postprocess_hawc2.py.
# IF ANYTHING GOES WRONG:
#   E.g., the job finishes, you refresh the folder but there is no stats file.
#   Read the lsf log files, especially the error file (paths above), to get more info.

# Print some basic information to LSF log file
date
START_TIME=`date +%s`
echo whoami
whoami
echo hostname -I
hostname -I

echo "Postprocessing HAWC2..."

# make the lsf_log folder if it doesn't exist
mkdir -p lsf_log/

# Load the Python module
module load python3/3.9.19

# Call Python on the post-processing script
python3 postprocess_hawc2.py


# Say goodbye
END_TIME=`date +%s`
echo done
echo "runtime: $(($END_TIME-$START_TIME)) seconds"