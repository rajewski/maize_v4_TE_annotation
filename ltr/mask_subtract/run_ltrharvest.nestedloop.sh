#!/bin/bash -l
#SBATCH -o ../../history/run_ltrharvest.nestedloop-%A.txt
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --mem=101G
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
if [ ! -e ${GENOMEBASE}.ltrharvest.contignames.gff3 ]; then
  echo "$(date +'%r'): Converting GFF3 from initial search to have proper contig names."
  python convert_ltrharvest_seq_gff_to_contignames.py ../${GENOMEBASE}.ltrharvest.gff3 > ${GENOMEBASE}.ltrharvest.contignames.gff3
  echo "$(date +'%r'): Done."
else
  echo "$(date +'%r'): Contig names from initial GFF3 have already been converted. Skipping to indexing."
fi

if [ ! -e ${GENOMEBASE}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 ]; then
  echo "$(date +'%r'): But first let me grep out the repeat regions from the GFF3, to bound the subtraction were about to do."
  grep --no-group-separator -B2 -A1 -P "LTR_retrotransposon\t" ${GENOMEBASE}.ltrharvest.contignames.gff3 | sed -n '1~2p' > ${GENOMEBASE}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3
  echo "$(date +'%r'): Done"
fi

if [ ! -e $GENOMEFASTA.2.fai ]; then
  echo "$(date +'%r'): Starting index creation with samtools."
  ## make an index for the r script and bedtools complement
  samtools faidx $GENOMEFASTA
  echo "$(date +'%r'): Index is done, now let's grab just the columns that bedtools will need."
  awk -v OFS='\t' '{print $1, $2}' $GENOMEFASTA.fai > $GENOMEFASTA.2.fai
  # mv $GENOMEFASTA.2.fai $GENOMEFASTA.fai #for some reason moving the file breaks this
else
  echo "$(date +'%r'): Indexing of genome is already complete. Skipping to complementing."
fi

if [ ! -e ${GENOMEBASE}.ltrharvest.contignames.NOTltrretrotransposon.gff3 ]; then
  echo "$(date +'%r'): Using bedtools to generate a list of where the LTRs are NOT located."
  bedtools complement -i ${GENOMEBASE}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 -g $GENOMEFASTA.2.fai > ${GENOMEBASE}.ltrharvest.contignames.NOTltrretrotransposon.gff3
  grep -Pv "scaffold\d*\t0\t0" $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.gff3 > $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.2.gff3 #fix the gff with 0 start and end lines
  awk '$2>$3{tmp = $3; $3=$2; $2=tmp;}; {print $1 "\t" $2 "\t" $3 }' $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.2.gff3 > $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.gff3 #flip columns that are backwards (http://www.metagenomics.wiki/tools/ubuntu-linux/awk)
  #mv $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.2.gff3 $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.gff3 #clean stuff up, not necessary because the awk command does this at the same time
  rm $GENOMEBASE.ltrharvest.contignames.NOTltrretrotransposon.2.gff3 #nut now this is necessary instead of the mv
  echo "$(date +'%r'): Done."
else
  echo  "$(date +'%r'): Complement already created. So let's make a FASTA file of those areas."
fi

if [ ! -e ${GENOMEBASE}.subtract1.fa ]; then
  ## generate a subtracted fasta
  echo "$(date +'%r'): Creating a subtracted FASTA for the next round of searching."
  bedtools getfasta -fi $GENOMEFASTA -bed ${GENOMEBASE}.ltrharvest.contignames.NOTltrretrotransposon.gff3 -fo ${GENOMEBASE}.subtract1.fa
  ## concatenate the entries by chromosome
  python collapse_chromosomes.py ${GENOMEBASE}.subtract1.fa > ${GENOMEBASE}.temp
  mv ${GENOMEBASE}.temp ${GENOMEBASE}.subtract1.fa
  echo "$(date +'%r'): Done."
else
  echo "$(date +'%r'): Subtracted FASTA already created. Let's move on to indexing it."
fi

## index this fasta
if [ ! -s $GENOMEBASE.subtract1.md5 ]; then
  echo "$(date +'%r'): Running Suffixerator to make genome index."
  $GENOMETOOLS suffixerator -db ${GENOMEBASE}.subtract1.fa -indexname ${GENOMEBASE}.subtract1 -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM
  echo "$(date +'%r'): Done."
else
  echo "$(date +'%r'): First round subtracted genome already indexed. Skipping to LTR Harvest."
fi

if [ ! -s subtract1/${GENOMEBASE}.subtract1.ltrharvest.gff3 ]; then
  echo "$(date +'%r'): Running LTR Harvest."
  mkdir -p subtract1
  $GENOMETOOLS ltrharvest -v -index ${GENOMEBASE}.subtract1 -gff3 subtract1/${GENOMEBASE}.subtract1.ltrharvest.gff3 -outinner subtract1/${GENOMEBASE}.subtract1.ltrharvest.outinner.fa -out subtract1/${GENOMEBASE}.subtract1.ltrharvest.fa > subtract1/${GENOMEBASE}.subtract1.ltrharvest.out
  echo "$(date +'%r'): Done."
  echo "$(date +'%r'): Sorting the GFF3 file."
  $GENOMETOOLS gff3 -sort subtract1/${GENOMEBASE}.subtract1.ltrharvest.gff3 > subtract1/${GENOMEBASE}.subtract1.ltrharvest.sorted.gff3
  echo "$(date +'%r'): Done."
else
  echo "$(date +'%r'): First round subtraction LTR Harvest output file found, so I'm skipping ahead to round 2."
fi

#Run the loop for 10 rounds
while [ $i -le 10 ]
do
  # Define some variables for the loop
  OLDINDEX=$i
  GENOME=${GENOMEBASE}.subtract${i}
  i=$(( $i + 1 ))
  NEWINDEX=$i
  NEWGENOME=${GENOMEBASE}.subtract${i}
  GENOMEFASTA=${GENOME}.fa
  NEWGENOMEFASTA=${NEWGENOME}.fa
  
  echo "$(date +'%r'): Running round ${NEWINDEX} subtraction."

  ### Run suffixerator to make a suffix array of the genome for genometools ###
  # switch genome tools back to their real contig names
  if [ ! -s subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.gff3 ]; then
    python convert_ltrharvest_seq_gff_to_contignames.py subtract${OLDINDEX}/${GENOME}.ltrharvest.gff3 > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.gff3
  fi
  
  ## only get the LTR_retrotransposon records to keep Rscript from repeating computation
  grep --no-group-separator -B2 -A1 -P "LTR_retrotransposon\t" subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.gff3 | sed -n '1~2p' > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3
  ## make an index for the r script and bedtools complement
  if [ ! -s ${GENOMEFASTA}.fa.fai ]; then
    echo "$(date +'%r'): Indexing ${GENOMEFASTA} for subtraction."
    samtools faidx ${GENOMEFASTA}
    awk -v OFS='\t' '{print $1, $2}' $GENOMEFASTA.fai > $GENOMEFASTA.2.fai
    echo "$(date +'%r'): Done."
  else
    echo "$(date +'%r'): ${GENOMEFASTA} already indexed. I'm skipping to subtracting from the gff3."
  fi
  
  if [ ! -s subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 ]; then
    echo "$(date +'%r'): Locating previously identified TEs in round $OLDINDEX gff3."
    ## find regions not covered by TEs
    bedtools complement -i subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 -g ${GENOME}.fa.2.fai > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3
    ## clean up the 0 lines in this file
    grep -Pv "\t0\t0" subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.1.gff3
    #awk '$2>$3{tmp = $3; $3=$2; $2=tmp;}; {print $1 "\t" $2 "\t" $3 }' subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.1.gff3 > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 #flip columns that are backwards (http://www.metagenomics.wiki/tools/ubuntu-linux/awk)
    #rm subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.1.gff3 # remove instead of mv
    mv subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.1.gff3 subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3
    echo "$(date +'%r'): Done."
  else
    echo "$(date +'%r'): Round ${OLDINDEX} subtraction TEs already located. Moving on to creating the ${NEWGENOME} fasta file."
  fi

  if [ ! -s $NEWGENOMEFASTA ]; then
    ## generate a subtracted fasta
    bedtools getfasta -fi $GENOMEFASTA -bed subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 -fo $NEWGENOMEFASTA
    ## concatenate the entries by chromosome
    python collapse_chromosomes.py $NEWGENOMEFASTA > ${NEWGENOMEFASTA}.temp
    mv ${NEWGENOMEFASTA}.temp $NEWGENOMEFASTA
    echo "$(date +'%r'): ${NEWGENOMEFASTA} successfully created and gaps from round ${OLDINDEX} TEs closed."
  else
    echo "$(date +'%r'): ${NEWGENOMEFASTA} already exists. Skipping to indexing this new fasta file."
  fi

  if [ ! -s ${NEWGENOMEFASTA}.md5 ]; then
    ### index this fasta
    echo "$(date +'%r'): Indexing ${NEWGENOMEFASTA}"
    $GENOMETOOLS suffixerator -db ${NEWGENOMEFASTA} -indexname ${NEWGENOME} -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM
    echo "$(date +'%r'): Done"
  else
    echo "$(date +'%r'): ${NEWGENOMEFASTA} already indexed. Skipping to running LTR harvest."
  fi

  ## Run LTR harvest
  if [ ! -s subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.sorted.gff3 ]; then
    echo "$(date +'%r'): Running LTR Harvest on ${NEWGENOMEFASTA}"
    mkdir -p subtract${NEWINDEX}
    # allow extra 1kb for each iteration, because we miss insertions that are not structural
    MAXLEN=$(($i * 1000 + 15000))    ### so 15kb default for the first round, plus the additional 1kb per round
    # all defaults except for maxdistltr 
    $GENOMETOOLS ltrharvest -v -index ${NEWGENOME} -gff3 subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.gff3 -outinner subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.outinner.fa -out subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.fa > subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.out
    $GENOMETOOLS gff3 -sort subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.gff3 > subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.sorted.gff3
    echo "$(date +'%r'): Done identifying TEs from ${NEWGENOMEFASTA}."
  else
    echo "$(date +'%r'): TEs from ${NEWGENOMEFASTA} already located. Moving on."
  fi
done

echo "$(date +'%r'): I'm done. Officially."