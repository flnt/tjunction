#!/bin/bash
#SBATCH --job-name=Tjunc
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 64
#SBATCH --cpus-per-task 1
#SBATCH --time 48:00:0

JVIEW="$1"
CASE_ARG="$2"
LENGTH="$3"
VELO="$4"
FORCE="$5"
LEVEL="$6"
MINLEVEL="$7"

OUTPUT_DIR="/scratch/coulombe/new_Tjunction_case${CASE_ARG}_l${LENGTH}_Ca${VELO}_f${FORCE}_LEVEL${LEVEL}"
mkdir -p "$OUTPUT_DIR"

cp Tjunction.c "$OUTPUT_DIR/"

cd "$OUTPUT_DIR" || exit

qcc -disable-dimensions -DDISPLAY=-$JVIEW -DTRACE=2 -source -D_MPI=1 Tjunction.c

module purge
module load gcc
module load openmpi
module load mesa
module load ffmpeg

mpicc -Wall -O2 -std=c99 -D_XOPEN_SOURCE=700 -DDISPLAY=-$JVIEW -DTRACE=2 _Tjunction.c -o Tjunction -I$BASILISK -L$BASILISK/gl -lglutils -lfb_tiny -L$BASILISK/wsServer -lws -lm

srun ./Tjunction $CASE_ARG $LENGTH $VELO $FORCE $LEVEL $MINLEVEL
