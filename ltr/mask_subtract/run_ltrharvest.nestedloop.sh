#SBATCH -o ../../history/run_ltrharvest.nestedloop-%A.txt
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=96G
#SBATCH --mail-user=araje002@ucr.edu
#SBATCH --mail-type=all
set -eu

# Load Software 
module load genometools/1.5.9
module load bedtools/2.25.0
module load samtools/1.6

# Set some variables to help
GENOMETOOLS=gt
MEMLIM=96GB
GENOMEBASE=NIOBT_r1.0
GENOMEFASTA=~/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.fasta
i=1

##### Set up first substraction becuase it needs some tweaks before we can make it routine
python convert_ltrharvest_seq_gff_to_contignames.py ../${GENOMEBASE}.ltrharvest.gff3 > ${GENOMEBASE}.ltrharvest.contignames.gff3
grep --no-group-separator -B2 -A1 "LTR_retrotransposon	" ${GENOMEBASE}.ltrharvest.contignames.gff3 | sed -n '1~2p' > ${GENOMEBASE}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3
## make an index for the r script and bedtools complement
samtools faidx $GENOMEFASTA
awk -v OFS='\t' '{print $1, $2}' $GENOMEFASTA.fai > $GENOMEFASTA.2.fai
mv $GENOMEFASTA.2.fai $GENOMEFASTA.fai
bedtools complement -i ${GENOMEBASE}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 -g $GENOMEFASTA.fai > ${GENOMEBASE}.ltrharvest.contignames.NOTltrretrotransposon.gff3

## generate a subtracted fasta
grep -Pv "scaffold\d*\t0\t0" $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.gff3 > $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.2.gff3 #fix the gff with 0 start and end lines
mv $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.2.gff3 $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.gff3 #clean stuff up
bedtools getfasta -fi $GENOMEFASTA -bed ${GENOMEBASE}.ltrharvest.contignames.NOTltrretrotransposon.gff3 -fo ${GENOMEBASE}.subtract1.fa

## concatenate the entries by chromosome
python collapse_chromosomes.py ${GENOMEBASE}.subtract1.fa > ${GENOMEBASE}.temp
mv ${GENOMEBASE}.temp ${GENOMEBASE}.subtract1.fa

## index this fasta
if [ ! -s $GENOMEBASE.substract1.md5 ]; then
  echo "Running Suffixerator to make genome index."
  $GENOMETOOLS suffixerator -db ${GENOMEBASE}.subtract1.fa -indexname ${GENOMEBASE}.subtract1 -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM
else
  echo "First round subtracted genome already indexed. Skipping to LTR Harvest."
fi

if [ ! -s subtract1/${GENOMEBASE}.subtract1.ltrharvest.gff3 ]; then
  echo "Running LTR Harvest."
  mkdir -p subtract1
  $GENOMETOOLS ltrharvest -index ${GENOMEBASE}.subtract1 -gff3 subtract1/${GENOMEBASE}.subtract1.ltrharvest.gff3 -outinner subtract1/${GENOMEBASE}.subtract1.ltrharvest.outinner.fa -out subtract1/${GENOMEBASE}.subtract1.ltrharvest.fa > subtract1/${GENOMEBASE}.subtract1.ltrharvest.out
  echo "Sorting the GFF3 file."
  $GENOMETOOLS gff3 -sort subtract1/${GENOMEBASE}.subtract1.ltrharvest.gff3 > subtract1/${GENOMEBASE}.subtract1.ltrharvest.sorted.gff3
else
  echo "First round subtraction LTR Harvest output file found, so I'm skipping ahead to round 2."
fi


### Let the while loop work for 5 rounds. My genome is fragmented, so the original 100 rounds seems excessive
#while [ grep -c ltr_retrotransposon ${GENOMEBASE}.hardmask${i}.ltrharvest.gff3 -gt 0 ]
while [ $i -le 5 ]
do
  # Define some variables for the loop
  OLDINDEX=$i
  GENOME=${GENOMEBASE}.subtract${i}
  i=$(( $i + 1 ))
  NEWINDEX=$i
  NEWGENOME=${GENOMEBASE}.subtract${i}
  GENOMEFASTA=${GENOME}.fa
  NEWGENOMEFASTA=${NEWGENOME}.fa
  
  echo "Running round ${NEWINDEX} subtraction."

  ### Run suffixerator to make a suffix array of the genome for genometools ###
  # switch genome tools back to their real contig names
  if [ ! -s subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.gff3 ]; then
    python convert_ltrharvest_seq_gff_to_contignames.py subtract${OLDINDEX}/${GENOME}.ltrharvest.gff3 > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.gff3
  fi
  
  ## only get the LTR_retrotransposon records to keep Rscript from repeating computation
  grep --no-group-separator -B2 -A1 "LTR_retrotransposon\t" subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.gff3 | sed -n '1~2p' > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3
  ## make an index for the r script and bedtools complement
  if [ ! -s ${GENOMEFASTA}.fa.fai ]; then
    echo "Indexing ${GENOMEFASTA} for subtraction."
    samtools faidx ${GENOMEFASTA}
    echo "Done."
  else
    echo "${GENOMEFASTA} already indexed. I'm skipping to subtracting from the gff3."
  fi
  
  if [ ! -s subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 ]; then
    echo "Locating previously identified TEs in round $OLDINDEX gff3."
    ## find regions not covered by TEs
    bedtools complement -i subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 -g ${GENOME}.fa.fai > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3
    ## clean up the extra columns in this file
    grep -Pv "scaffold\d*\t0\t0" subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.1.gff3
    mv subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.1.gff3 subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3
    echo "Done."
  else
    echo "Round ${OLDINDEX} subtraction TEs already located. Moving on to creating the ${NEWGENOME} fasta file."
  fi

  if [ ! -s $NEWGENOMEFASTA ]; then
    ## generate a subtracted fasta
    bedtools getfasta -fi $GENOMEFASTA -bed subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 -fo $NEWGENOMEFASTA
    ## concatenate the entries by chromosome
    python collapse_chromosomes.py $NEWGENOMEFASTA > ${NEWGENOMEFASTA}.temp
    mv ${NEWGENOMEFASTA}.temp $NEWGENOMEFASTA
    echo "${NEWGENOMEFASTA} successfully created and gaps from round ${OLDINDEX} TEs closed."
  else
    echo "${NEWGENOMEFASTA} already exists. Skipping to indexing this new fasta file."
  fi

  if [ ! -s ${NEWGENOMEFASTA}.md5 ]; then
    ### index this fasta
    echo "Indexing ${NEWGENOMEFASTA}"
    $GENOMETOOLS suffixerator -db ${NEWGENOMEFASTA} -indexname ${NEWGENOME} -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM
    echo "Done"
  else
    echo "${NEWGENOMEFASTA} already indexed. Skipping to running LTR harvest."
  fi

  ## Run LTR harvest
  if [ ! -s subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.sorted.gff3 ]; then
    echo "Running LTR Harvest on ${NEWGENOMEFASTA}"
    mkdir -p subtract${NEWINDEX}
    # allow extra 1kb for each iteration, because we miss insertions that are not structural
    MAXLEN=$(($i * 1000 + 15000))    ### so 15kb default for the first round, plus the additional 1kb per round
    # all defaults except for maxdistltr 
    $GENOMETOOLS ltrharvest -index ${NEWGENOME} -gff3 subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.gff3  -maxdistltr $MAXLEN -outinner subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.outinner.fa -out subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.fa > subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.out
    $GENOMETOOLS gff3 -sort subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.gff3 > subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.sorted.gff3
    echo "Done identifying TEs from ${NEWGENOMEFASTA}."
  else
    echo "TEs from ${NEWGENOMEFASTA} already located. Moving on."
done



# 10 Nov 18, this is probably not necessary below:

## can then run all the ltrdigest in array form on all the gffs
###################
## run ltrdigest ##
###################

#mkdir -p hardmask5/ltrdigest

#$GENOMETOOLS -j 16 ltrdigest -outfileprefix hardmask5/ltrdigest/$GENOME.ltrdigest -trnas eukaryotic-tRNAs.fa -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -- hardmask5/$GENOME.ltrharvest.sorted.gff3 $GENOME > hardmask5/$GENOME.ltrdigest.gff3

