#!/usr/bin/env bash

## v3.0 of STARsolo wrappers is set up to guess the chemistry automatically
## newest version of the script uses at least STAR v2.7.9a with EM multimapper processing in STARsolo (it's on by default; turn this off deleting --soloMultiMappers EM)
## velocyto extra processing also became unnecessary
echo "Aligning and counting RNA-velocity matrixes for sample ${snakemake_wildcards[run]} with STARsolo" 2> "${snakemake_log[0]}"

TAG=${snakemake_params[sample_run_name]}
if [[ $TAG == "" ]]
then
  echo "Usage: ./starsolo_10x_auto.sh <sample_tag>" 2>> "${snakemake_log[0]}"
  echo "(make sure you set the correct REF, FQDIR, and SORTEDBAM/NOBAM variables)" 2>> "${snakemake_log[0]}"
  exit 1
fi

CWD=`pwd`
CPUS=${snakemake[threads]}                                             ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.
REF="${snakemake_input[reference]}"                                    ## choose the appropriate reference
WL="${snakemake_input[barcodes]}"                                      ## directory with all barcode whitelists
FQDIR="$CWD/${snakemake_input[sample]}"
OUTDIR="$CWD/${snakemake_output[sample]}"

R1="$CWD/${snakemake_input[r1]}"
R2="$CWD/${snakemake_input[r2]}"

## choose one of the two otions, depending on whether you need a BAM file
#BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --limitBAMsortRAM 120000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB GX GN"
BAM="--outSAMtype None"

###################################################################### DONT CHANGE OPTIONS BELOW THIS LINE ##############################################################################################

if [[ `which samtools` == "" || `which seqtk` == "" || `which STAR` == "" ]]
then
  echo "ERROR: Please make sure you have STAR (v2.7.9a or above), samtools, and seqtk installed and available in PATH!" 2>> "${snakemake_log[0]}"
  exit 1
fi

## let's see if the files are archived or not. Gzip is the only one we test for, but bgzip archives should work too since they are gzip-compatible.
GZIP=""
BC=""
NBC1=""
NBC2=""
NBC3=""
NBCA=""
R1LEN=""
R2LEN=""
R1DIS=""

mkdir -p $OUTDIR &&  cd $OUTDIR

## randomly subsample 200k reads - let's hope there are at least this many (there should be):
seqtk sample -s100 $R1 200000 > test.R1.fastq &
seqtk sample -s100 $R2 200000 > test.R2.fastq &
wait

## see if the original fastq files are archived:
if [[ `find $FQDIR/* -maxdepth 1 | grep $TAG | grep "\.gz$"` != "" ]]
then
  GZIP="--readFilesCommand zcat"
fi

NBC1=`cat test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/737K-april-2014_rc.txt | wc -l`
NBC2=`cat test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/737K-august-2016.txt | wc -l`
NBC3=`cat test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/3M-february-2018.txt | wc -l`
NBCA=`cat test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/737K-arc-v1.txt | wc -l`
R1LEN=`cat test.R1.fastq | awk 'NR%4==2' | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}'`
R2LEN=`cat test.R2.fastq | awk 'NR%4==2' | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}'`
R1DIS=`cat test.R1.fastq | awk 'NR%4==2' | awk '{print length($0)}' | sort | uniq -c | wc -l`

## elucidate the right barcode whitelist to use. Grepping out N saves us some trouble. Note the special list for multiome experiments (737K-arc-v1.txt):
if (( $NBC2 > 100000 ))
then
  BC=$WL/737K-august-2016.txt
elif (( $NBC3 > 100000 ))
then
  BC=$WL/3M-february-2018.txt
elif (( $NBCA > 100000 ))
then
  BC=$WL/737K-arc-v1.txt
elif (( $NBC1 > 100000 ))
then
  BC=$WL/737K-april-2014_rc.txt
else
  >&2 echo "ERROR: No whitelist has matched a random selection of 200,000 barcodes!"
  exit 1
fi

## check read lengths, fail if something funky is going on:
PAIRED=False
UMILEN=""
CBLEN=""
if (( $R1DIS > 1 && $R1LEN <= 30 ))
then
  echo "ERROR: Read 1 (barcode) has varying length; possibly someone thought it's a good idea to quality-trim it. Please check the fastq files." 2>> "${snakemake_log[0]}"
  exit 1
elif (( $R1LEN < 24 ))
then
  echo "ERROR: Read 1 (barcode) is less than 24 bp in length. Please check the fastq files." 2>> "${snakemake_log[0]}"
  exit 1
elif (( $R2LEN < 40 ))
then
  echo "ERROR: Read 2 (biological read) is less than 40 bp in length. Please check the fastq files." 2>> "${snakemake_log[0]}"
  exit 1
fi

## assign the necessary variables for barcode/UMI length/paired-end processing:
if (( $R1LEN > 50 ))
then
  PAIRED=True
  UMILEN=10
  CBLEN=16
elif (( $NBC1 > 100000 ))
then
  CBLEN=14
  UMILEN=$((R1LEN-14))
else
  CBLEN=16
  UMILEN=$((R1LEN-16))
fi

## finally, see if you have 5' or 3' experiment. I don't know and easier way:
STRAND=Forward

STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn test.R2.fastq test.R1.fastq --runDirPerm All_RWX --outSAMtype None --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand Forward --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloFeatures Gene GeneFull --soloOutFileNames test_strand/ features.tsv barcodes.tsv matrix.mtx &> /dev/null 2>> "${snakemake_log[0]}"

## the following is needed in case of bad samples: when a low fraction of reads come from mRNA, experiment will look falsely reverse-stranded
UNIQFRQ=`grep "Reads Mapped to Genome: Unique," test_strand/GeneFull/Summary.csv | awk -F "," '{print $2}'`
GENEPCT=`grep "Reads Mapped to GeneFull: Unique GeneFull" test_strand/GeneFull/Summary.csv | awk -F "," -v v=$UNIQFRQ '{printf "%d\n",$2*100/v}'`

if (( $GENEPCT < 20 ))
then
  STRAND=Reverse
fi

## finally, if paired-end experiment turned out to be 3' (yes, they do exist!), process it as single-end:
if [[ $STRAND == "Forward" && $PAIRED == "True" ]]
then
  PAIRED=False
  if [[ $BC == "$WL/3M-february-2018.txt" ]]
  then
    UMILEN=12
  fi
fi

echo "Done setting up the STARsolo run; here are final processing options:" 2>> "${snakemake_log[0]}"
echo "=============================================================================" 2>> "${snakemake_log[0]}"
echo "Sample: $TAG" 2>> "${snakemake_log[0]}"
echo "Paired-end mode: $PAIRED" 2>> "${snakemake_log[0]}"
echo "Strand (Forward = 3', Reverse = 5'): $STRAND, %reads same strand as gene: $GENEPCT" 2>> "${snakemake_log[0]}"
echo "CB whitelist: $BC, matches out of 200,000: $NBC3 (v3), $NBC2 (v2), $NBC1 (v1), $NBCA (multiome) " 2>> "${snakemake_log[0]}"
echo "CB length: $CBLEN" 2>> "${snakemake_log[0]}"
echo "UMI length: $UMILEN" 2>> "${snakemake_log[0]}"
echo "GZIP: $GZIP" 2>> "${snakemake_log[0]}"
echo "-----------------------------------------------------------------------------" 2>> "${snakemake_log[0]}"
echo "Read 1 files: $R1" 2>> "${snakemake_log[0]}"
echo "-----------------------------------------------------------------------------" 2>> "${snakemake_log[0]}"
echo "Read 2 files: $R2" 2>> "${snakemake_log[0]}"
echo "-----------------------------------------------------------------------------" 2>> "${snakemake_log[0]}"

if [[ $PAIRED == "True" ]]
then
  STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R1 $R2 --runDirPerm All_RWX $GZIP $BAM --soloBarcodeMate 1 --clip5pNbases 39 0 --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloCBstart 1 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand Forward --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR --outFilterScoreMin 30 --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx --soloMultiMappers EM --outReadsUnmapped Fastx  2>> "${snakemake_log[0]}"
else
  STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R2 $R1 --runDirPerm All_RWX $GZIP $BAM --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand $STRAND --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx --soloMultiMappers EM --outReadsUnmapped Fastx  2>> "${snakemake_log[0]}"
fi

## index the BAM file
if [[ -s Aligned.sortedByCoord.out.bam ]]
then
  samtools index -@$CPUS Aligned.sortedByCoord.out.bam
fi

## finally, let's gzip all outputs
gzip Unmapped.out.mate1 &
gzip Unmapped.out.mate2 &

cd output
for i in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered Velocyto/raw Velocyto/filtered
do
  cd $i; for j in *; do gzip $j & done
  cd ../../
done

wait
echo "ALL DONE!" 2>> "${snakemake_log[0]}"
