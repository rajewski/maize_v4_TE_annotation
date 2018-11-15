#!/bin/bash -l
#SBATCH -o ../../history/run_ltrdigest.array-%A-%a.txt
#SBATCH --ntasks=30
#SBATCH --nodes=1
#SBATCH --time=2-20:00:00
#SBATCH --mem=101G
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=all
set -eu

# Set some variables
module load genometools/1.5.9
module load hmmer/3.1b2
GENOMETOOLS=gt
GENOMEBASE=NIOBT_r1.0
i=$SLURM_ARRAY_TASK_ID
GENOME=${GENOMEBASE}.subtract${i}
GENOMEFASTA=${GENOME}.fa
MEMLIM=96GB

## run ltrdigest
mkdir -p subtract${i}/ltrdigest
$GENOMETOOLS -j $SLURM_NTASKS ltrdigest -outfileprefix subtract${i}/ltrdigest/$GENOME.ltrdigest -trnas ../eukaryotic-tRNAs/eukaryotic-tRNAs.fa -hmms ../gydb_hmms/GyDB_collection/profiles/*.hmm -- subtract${i}/$GENOME.ltrharvest.sorted.gff3 $GENOME > subtract${i}/$GENOME.ltrdigest.gff3


### fix the .des files that look funny
#gt encseq encode -des yes -dna yes B73.Mhap2.quiver.subtract13.fa
#$GENOMETOOLS encseq encode -des yes -ssp no -sds no -md5 no -dna yes -indexname $GENOME $GENOMEFASTA

#$GENOMETOOLS suffixerator -db ${GENOMEFASTA} -indexname ${GENOME}.try -des -ssp no -sds no -md5 no -dna -memlimit $MEMLIM

