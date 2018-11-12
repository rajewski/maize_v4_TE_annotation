#!/bin/bash -login
#SBATCH -o ../../history/target_mutator-%A-%a.txt
#SBATCH --ntasks=4
#SBATCH --mem=6G
#SBATCH --time=00:45:00
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=all
#SBATCH --array=0-137
set -eu

### Delcare some variables and load software
GENOME=~/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.fasta
module load TARGeT/2.10
module unload ncbi-blast
module load ncbi-blast/2.2.26
export PATH=$PATH:$HOME/software/mTEA/lib/blogo/
export PERL5LIB=$PERL5LIB:$HOME/software/mTEA/lib/blogo/mTEA/lib/blogo/
#export PATH=$PATH:$HOME/software/ImageMagick-7.0.8-14/ #maybe not necessary

### set up a file name for each individual TE fasta
FILES=($(ls -1 $HOME/TE_Annotation/tir/target_mTEA/tir_references/DTM_*))
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
    target.py -q $FILENAME -t nucl -o $DNATE -i s -P $SLURM_NTASKS -S PHI -b_a 10000 -b_d 10000 -p_n 10000 -p_f 200 -p_M 0.2 $GENOME ${DNATE}_target
    echo "$(date +'%r'): Done."
  else
    echo "$(date +'%r'): TARGeT index database of genome already created."
    echo "$(date +'%r'): Running TARGeT for ${DNATE}."
    target.py -q $FILENAME -t nucl -o $DNATE -i s -P $SLURM_NTASKS -S PHI -DB -b_a 10000 -b_d 10000 -p_n 10000 -p_f 200 -p_M 0.2 $GENOME ${DNATE}_target
    echo "$(date +'%r'): Done."
  fi
  ### Convert names of flanking fasta file
  if [ ! -f ${DNATE}.flank.fa ]; then
    python convert_target_toTIRID.py ${DNATE}/*/${DNATE}.flank > ${DNATE}.flank.fa
  else
    python convert_target_toTIRID.py ${DNATE}/*/${DNATE}.flank_adj > ${DNATE}.flank.fa
  fi
  if [ ! -f ${DNATE}.tir.fa ]; then
    ## Mutator : -c NNNNNNNNN -t NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN -d 8 -s 1
    ## (only 40 bp TIR because Lisch and Jiang talk about non-TIR MULEs with TIR <50bp. also worried 200bp will not blast well)
    echo "$(date +'%r'): Running mTEA on ${DNATE}."
    perl id_TIR_in_FASTA.mcs.pl -o ${DNATE}.tir.fa -i ${DNATE}.flank.fa -c NNNNNNNNN -t NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN -d 8 -s 1
    echo "$(date +'%r'): Done."
  else
    echo "$(date +'%r'): mTEA has already been run for ${DNATE}. Exiting."
  fi
  if [ ! -s ${DNATE}/*/${DNATE}.flank ]; then
    echo "$(date +'%r'): No matches to ${DNATE} found. Deleting empty files."
    rm -rf $DNATE*
  fi
else
  echo "$(date +'%r'): TARGeT has already been run on ${DNATE}. Exiting."
fi