#!/bin/bash -login
#SBATCH -o ../history/sine-%A.txt
#SBATCH -J sine_finder
#SBATCH --ntasks=16
#SBATCH --mem=96G
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=all
set -eu

GENOME=NIOBT_r1.0
GENOMEFASTA=~/shared/Nobtusifolia/Genome_Files/${GENOME}.fasta

# Load modules/software
export PATH=$PATH:$HOME/software/silix-1.2.9/bin
SILIX=silix
module load vsearch/1.10.2
VSEARCH=vsearch

# Download SINE-FINDER?
if [ ! -x sine_finder.py ]; then
  echo "Downloading SINE-FINDER from plantcell.org."
  # get the SINE-Finder program from Wenke et al. 2011
  wget http://www.plantcell.org/content/suppl/2011/08/29/tpc.111.088682.DC1/Supplemental_Data_Set_1-sine_finder.txt 
  # rename ot python script
  mv Supplemental_Data_Set_1-sine_finder.txt sine_finder.py
  # make executable
  chmod 774 sine_finder.py
else
  echo "SINE-FINDER already downloaded, so I'm going to run it now."
fi

# Run sinefinder
#### I haven't been able to get sine_finder to work with reverse sequences, as it seems to report TSDs wrong on the reverse strand.
####   so I'm only reporting on the forward strand.
### -f both : outputs csv and fasta
if [ ! -f ${GENOME}-matches.fasta ]; then
  echo "Running SINE-FINDER."
  python sine_finder.py -T chunkwise -V1 -f both ${GENOMEFASTA}
  mv ~/shared/Nobtusifolia/Genome_Files/${GENOME}-matches.fasta .
  mv ~/shared/Nobtusifolia/Genome_Files/${GENOME}-matches.csv .
  echo "Done."
  if [ ! -f ${GENOME}-matches.noTSD.fa }; then
    echo "Removing TSDs from the results."
    python remove_tsd_sinefinder.py ${GENOME}-matches.fasta ${GENOME}-matches.noTSD.fa
    echo "Done"
  else
    #this case should never happen, but whatever
    echo "TSDs have already been removed from the results. Skipping to clustering."
  fi
else
  echo "SINE-FINDER has already been run."
  if [ ! -f ${GENOME}-matches.noTSD.fa }; then
    echo "But the TSDs haven't been removed from the results. I'll do that now."
    python remove_tsd_sinefinder.py ${GENOME}-matches.fasta ${GENOME}-matches.noTSD.fa
    echo "Done"
  else
    echo "And TSDs have also already been removed from the results. Skipping to clustering of results."
  fi
fi


#### vsearch to identify homology, silix to cluster
$VSEARCH -allpairs_global ${GENOME}-matches.noTSD.fa -blast6out ${GENOME}-matches.noTSD.allvall.8080.out -id 0.8 -query_cov 0.8 -target_cov 0.8 --threads $SLURM_NTASKS

# single linkage cluster those that are 80% identical to each other.
$SILIX ${GENOME}-matches.noTSD.fa ${GENOME}-matches.noTSD.allvall.8080.out -f SINE -i 0.8 -r 0.8 > ${GENOME}-matches.noTSD.8080.fnodes

#This clustering of an external DB no longer works because MTEC is down.
# cluster my families into MTEC TE families
#wget http://maizetedb.org/~maize/TE_12-Feb-2015_15-35.fa
#$VSEARCH --usearch_global TE_12-Feb-2015_15-35.fa -db ${GENOME}-matches.noTSD.fa -id 0.8 -query_cov 0.8 -target_cov 0.8 -blast6out ${GENOME}-matches.noTSD.TEDB8080.out -strand both -top_hits_only --threads $SLURM_NTASKS

### cluster into families and output final gff with this R script
# Rscript generate_gff_SINE.R $GENOME #I am going to run this separately the first time, to make it work without the MTEC clusters


