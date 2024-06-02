#!/bin/bash
#SBATCH --job-name=inlet
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 64
#SBATCH --cpus-per-task 1
#SBATCH --time 48:00:0

JVIEW="$1"
Re="$2"
Ca="$3"
LI="$4"
L="$5"
LEVEL="$6"
MINLEVEL="$7"

OUTPUT_DIR="/scratch/coulombe/new_inlet_Re${Re}_Ca${Ca}_LI${LI}_L${L}_LEVEL${LEVEL}"
mkdir -p "$OUTPUT_DIR"

cp inlet.c "$OUTPUT_DIR/"

cd "$OUTPUT_DIR" || exit

qcc -disable-dimensions -DDISPLAY=-$JVIEW -source -D_MPI=1 inlet.c

module purge
module load gcc
module load openmpi
module load mesa
module load ffmpeg

mpicc -Wall -O2 -std=c99 -D_XOPEN_SOURCE=700 -DDISPLAY=-$JVIEW _inlet.c -o inlet -I$BASILISK -L$BASILISK/gl -lglutils -lfb_tiny -L$BASILISK/wsServer -lws -lm

srun ./inlet $Re $Ca $LI $L $LEVEL $MINLEVEL
