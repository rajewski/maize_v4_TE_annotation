#!/bin/bash -login
#SBATCH -o ../history/target_cacta-%A.txt
#SBATCH --ntasks=16
#SBATCH --mem=24G
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=all
set -eu

### Delcare some variables and load software
GENOME=~/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.fasta
module load TARGeT/2.10
module unload ncbi-blast
module load ncbi-blast/2.2.26
export PATH=$PATH:$HOME/software/mTEA/lib/blogo/:$HOMEsoftware/ImageMagick-7.0.8-14
export PERL5LIB=$PERL5LIB:$HOME/software/mTEA/lib/blogo/mTEA/lib/blogo/

### set up a file name for each individual TE fasta
FILES=($(ls -1 $HOME/TE_Annotation/tir/target_mTEA/tir_references/DTC_*))
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}
echo $FILENAME

### RUN TARGeT
FILE=$(basename "$FILENAME")
echo $FILE
DNATE="${FILE%.*}"
echo $DNATE

if [ ! -f ${DNATE}.tir.fa.tab ]; then
  echo "$(date +'%r'): Running TARGeT on ${DNATE}."
  mkdir -p $DNATE
  if [ ! -s ${GENOME}.index ]; then
    ### db generation for the first time to index the genome
    echo "$(date +'%r'): Generating TARGeT index database of the genome and running TARGeT for ${DNATE}."
    target.py -q $FILENAME -t nucl -o $DNATE -i s -P 1 -S PHI -b_a 10000 -b_d 10000 -p_n 10000 -p_f 200 -p_M 0.3 $GENOME ${DNATE}_target -T 1
    echo "$(date +'%r'): Done."
  else
    echo "$(date +'%r'): TARGeT index database of genome already created."
    echo "$(date +'%r'): Running TARGeT for ${DNATE}."
    target.py -q $FILENAME -t nucl -o $DNATE -i s -P 1 -S PHI -DB -b_a 10000 -b_d 10000 -p_n 10000 -p_f 200 -p_M 0.3 $GENOME ${DNATE}_target
    echo "$(date +'%r'): Done."
  fi
  ### Convert names of flanking fasta file
  if [ ! -f ${DNATE}.flank_adj ]; then
    python convert_target_toTIRID.py ${DNATE}/*/${DNATE}.flank > ${DNATE}.flank.fa
  else
    python convert_target_toTIRID.py ${DNATE}/*/${DNATE}.flank_adj > ${DNATE}.flank.fa
  fi
  if [ ! -f ${DNATE}.tir.fa ]; then
    ### RUN mTEA for CACTA
    echo "$(date +'%r'): Running mTEA on ${DNATE}."
    perl id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c NNN -t CACTNNNNNNNNN -d 3 -s 0
    echo "$(date +'%r'): Done."
    #rm -rf $DNATE
  else
    echo "$(date +'%r'): mTEA has already been run for ${DNATE}. Exiting."
fi

### so to submit, grab files and export them
####### DO ALL OF THIS IN THE SHELL YOU'RE ABOUT TO SUBMIT IN
### export FILES=($(ls -1 DTA*))
### NUMFILES=${#FILES[@]}
### ARRAYNUM=$(($NUMFILES - 1)) ## because slurm array is zero based

## if [ $ARRAYNUM -ge 0 ]; then
## sbatch --array=0-$ARRAYNUM tir_search_genome.sh
## fi
