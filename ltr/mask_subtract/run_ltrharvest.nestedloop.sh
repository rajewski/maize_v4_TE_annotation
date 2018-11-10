#SBATCH -o ../../history/ltrharvest.nestedloop-%A.txt
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
GENOME=NIOBT_r1.0
GENOMEFASTA=~/shared/Nobtusifolia/Genome_Files/NIOBT_r1.0.fasta
i=1

##### Set up first substraction becuase it needs some tweaks before we can make it routine
python convert_ltrharvest_seq_gff_to_contignames.py ../${GENOME}.ltrharvest.gff3 > ${GENOME}.ltrharvest.contignames.gff3
grep --no-group-separator -B2 -A1 "LTR_retrotransposon	" ${GENOME}.ltrharvest.contignames.gff3 | sed -n '1~2p' > ${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3
## make an index for the r script and bedtools complement
samtools faidx $GENOMEFASTA
awk -v OFS='\t' '{print $1, $2}' $GENOMEFASTA.fai > $GENOMEFASTA.2.fai
mv $GENOMEFASTA.2.fai $GENOMEFASTA.fai
bedtools complement -i ${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 -g $GENOMEFASTA.fai > ${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3

## generate a subtracted fasta
grep -Pv "scaffold\d*\t0\t0" $GENOME.ltrharvest.contignames.NOTltrretrotransposon.gff3 > $GENOME.ltrharvest.contignames.NOTltrretrotransposon.2.gff3 #fix the gff with 0 start and end lines
mv $GENOME.ltrharvest.contignames.NOTltrretrotransposon.2.gff3 $GENOME.ltrharvest.contignames.NOTltrretrotransposon.gff3 #clean stuff up
bedtools getfasta -fi $GENOMEFASTA -bed ${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 -fo ${GENOME}.subtract1.fa

## concatenate the entries by chromosome
python collapse_chromosomes.py ${GENOME}.subtract1.fa > ${GENOME}.temp
mv ${GENOME}.temp ${GENOME}.subtract1.fa

## index this fasta
if [ ! -e $GENOME.substract1.md5 ]; then
  echo "Running Suffixerator to make genome index."
  $GENOMETOOLS suffixerator -db ${GENOME}.subtract1.fa -indexname ${GENOME}.subtract1 -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM
else
  echo "First round subtracted genome already indexed. Skipping to LTR Harvest."
fi

if [ ! -s subtract1/${GENOME}.subtract1.ltrharvest.out ]; then
  echo "Running LTR Harvest."
  mkdir -p subtract1
  $GENOMETOOLS ltrharvest -index ${GENOME}.subtract1 -gff3 subtract1/${GENOME}.subtract1.ltrharvest.gff3 -outinner subtract1/${GENOME}.subtract1.ltrharvest.outinner.fa -out subtract1/${GENOME}.subtract1.ltrharvest.fa > subtract1/${GENOME}.subtract1.ltrharvest.out
  echo "Sorting the GFF3 file."
  $GENOMETOOLS gff3 -sort subtract1/${GENOME}.subtract1.ltrharvest.gff3 > subtract1/${GENOME}.subtract1.ltrharvest.sorted.gff3
else
  echo "First round subtraction LTR Harvest output file found, so I'm skipping ahead to round 2."
fi


### go until there are no more LTR TEs in this nested form 
#while [ grep -c ltr_retrotransposon ${GENOMEBASE}.hardmask${i}.ltrharvest.gff3 -gt 0 ]
while [ $i -le 100 ]
do

# name the file stem based on suffixator index

OLDINDEX=$i
GENOME=${GENOMEBASE}.subtract${i}
i=$(( $i + 1 ))
NEWINDEX=$i
NEWGENOME=${GENOMEBASE}.subtract${i}

GENOMEFASTA=${GENOME}.fa
NEWGENOMEFASTA=${NEWGENOME}.fa

MEMLIM=96GB
CPU=16

###########################################################################
## Run suffixerator to make a suffix array of the genome for genometools ##
###########################################################################


## switch genome tools back to their real contig names
python convert_ltrharvest_seq_gff_to_contignames.py subtract${OLDINDEX}/${GENOME}.ltrharvest.gff3 > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.gff3
## only get the LTR_retrotransposon records to keep Rscript from repeating computation
grep --no-group-separator -B2 -A1 "LTR_retrotransposon	" subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.gff3 | sed -n '1~2p' > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3
## make an index for the r script and bedtools complement
samtools faidx ${GENOMEFASTA}

### run the rscript: read in gff, read in RDS of genomelist, update genomelist with updatePos, write gff with changed positions


## find regions not covered by TEs
bedtools complement -i subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.tsd.ltrretrotransposon.gff3 -g ${GENOME}.fa.fai > subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3

## generate a subtracted fasta
bedtools getfasta -fi $GENOMEFASTA -bed subtract${OLDINDEX}/${GENOME}.ltrharvest.contignames.NOTltrretrotransposon.gff3 -fo $NEWGENOMEFASTA

## concatenate the entries by chromosome
python collapse_chromosomes.py $NEWGENOMEFASTA > ${NEWGENOMEFASTA}.temp
mv ${NEWGENOMEFASTA}.temp $NEWGENOMEFASTA

### index this fasta
$GENOMETOOLS suffixerator -db ${NEWGENOMEFASTA} -indexname ${NEWGENOME} -tis -suf -lcp -des -ssp -sds -dna -memlimit $MEMLIM


#####################
## Run LTR harvest ##
#####################

mkdir -p subtract${NEWINDEX}
### allow extra 1kb for each iteration, because we miss insertions that are not structural
MAXLEN=$(($i * 1000 + 20000))    ### so 20kb for the first hardmask, plus the additional 1kb per round
## all defaults except for maxdistltr (default 15000)
$GENOMETOOLS ltrharvest -index ${NEWGENOME} -gff3 subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.gff3 -motif tgca -minlenltr 100 -maxlenltr 7000 -mindistltr 1000 -maxdistltr $MAXLEN -similar 85 -motifmis 1 -mintsd 5 -xdrop 5 -overlaps best -longoutput -outinner subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.outinner.fa -out subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.fa > subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.out

$GENOMETOOLS gff3 -sort subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.gff3 > subtract${NEWINDEX}/${NEWGENOME}.ltrharvest.sorted.gff3


done


## can then run all the ltrdigest in array form on all the gffs
###################
## run ltrdigest ##
###################

#mkdir -p hardmask5/ltrdigest

#$GENOMETOOLS -j 16 ltrdigest -outfileprefix hardmask5/ltrdigest/$GENOME.ltrdigest -trnas eukaryotic-tRNAs.fa -hmms gydb_hmms/GyDB_collection/profiles/*.hmm -- hardmask5/$GENOME.ltrharvest.sorted.gff3 $GENOME > hardmask5/$GENOME.ltrdigest.gff3

