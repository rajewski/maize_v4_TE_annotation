#!/bin/bash -login
#SBATCH -o ../history/ltrharvest-%A.txt
#SBATCH --ntasks=5
#SBATCH --nodes=1
#SBATCH --mem=96G
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=all
set -eu

# Load Software 
module load genometools/1.5.9
module load hmmer/3.1b2

# Set some variables to help
GENOMETOOLS=gt
MEMLIM=96GB
GENOME=NIOBT_r1.0
GENOMEFASTA=~/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.fasta

###########################################################################
## Run suffixerator to make a suffix array of the genome for genometools ##
###########################################################################

if [ ! -d index ]; then
  echo "Index directory not found, so I'll create it."
  mkdir index
  echo "Running Suffixerator to make genome index."
  $GENOMETOOLS suffixerator -db $GENOMEFASTA -indexname index/$GENOME -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM
  #make a symlink to the genomefasta so that and the index are in the same folder
  echo "Symlinking a copy of the genome fasta into the index folder."
  ln -s $GENOMEFASTA index/$GENOME.fasta
else
  if [ ! -e index/$GENOME.md5 ]; then
    echo "Running Suffixerator to make genome index."
    $GENOMETOOLS suffixerator -db $GENOMEFASTA -indexname index/$GENOME -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM
    #make a symlink to the genomefasta so that and the index are in the same folder
    echo "Symlinking a copy of the genome fasta into the index folder."
    ln -s $GENOMEFASTA index/$GENOME.fasta
  else
    echo "Index already made. Skipping to LTR Harvest"
  fi
fi

#####################
## Run LTR harvest ##
#####################

if [ ! -s ${GENOME}.ltrharvest.out ]; then
  echo "Running LTR Harvest."
  mkdir -p outinner
  $GENOMETOOLS ltrharvest -index index/$GENOME -gff3 $GENOME.ltrharvest.gff3 -outinner outinner/${GENOME}.ltrharvest.outinner.fa -out ${GENOME}.ltrharvest.fa > ${GENOME}.ltrharvest.out
  echo "Sorting the GFF3 file."
  $GENOMETOOLS gff3 -sort $GENOME.ltrharvest.gff3 > $GENOME.ltrharvest.sorted.gff3
else
  echo "LTR Harvest output file found, so I'm skipping to LTR Digest."
fi

###################
## run ltrdigest ##
###################

if [ ! -s $GENOME.ltrdigest.gff3 ]; then
  echo "Running LTR Digest."
  mkdir -p ltrdigest
  $GENOMETOOLS -j $SLURM_NTASKS ltrdigest -outfileprefix ltrdigest/$GENOME.ltrdigest -trnas eukaryotic-tRNAs/eukaryotic-tRNAs.fa -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -- $GENOME.ltrharvest.sorted.gff3 index/$GENOME > $GENOME.ltrdigest.gff3
else
  echo "LTR Digest output file found, so I'm done."
fi




