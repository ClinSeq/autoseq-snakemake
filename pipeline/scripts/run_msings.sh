#!/bin/bash

# The call of the script
echo "$0 $@"
echo

# Check if enough number of arguments
if [ $# -ne 11 ]; then
echo "Invalid number of options/arguments.
Usage: $0 [options] in.bam
Options (all required):
-b, --bed FILE          msi bed file, specifying the microsatellite regions for the utilized panel
-f, --fasta FILE        fasta reference genome
-i, --intervals FILE    msi intervals file customized for the utilized panel, for internal program use
-n, --baseline FILE     msi baseline file, based on msi negative samples with the utilized panel
-o, --outdir DIR        output directory" 1>&2
exit 1
fi

# parse the options
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -b|--bed)  # msi bed file
    BEDFILE="$2"
    shift # past argument
    shift # past value
    ;;
    -f|--fasta)  # fasta reference genome
    REF_GENOME="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--intervals)  # msi intervals file
    INTERVALS_FILE="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--baseline)  # msi baseline file
    MSI_BASELINE="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--outdir)  # output directory
    OUTDIR="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option / positional argument
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

# bam file to analyze, given as positional argument 1
BAM=$1

# activate the msings virtual environment
source $MSINGSENV/bin/activate  # NB! Need to set up the msings variable in alascca-dotfiles
VARSCAN=$MSINGSENV/bin/VarScan.v2.3.7.jar  # NB! Need to set up the msings variable in alascca-dotfiles


#"multiplier" is the number of standard deviations from the baseline that is required to call instability
multiplier=2.0 
#"msi_min_threshold" is the maximum fraction of unstable sites allowed to call a specimen MSI negative     
msi_min_threshold=0.1
#"msi_max_threshold" is the minimum fraction of unstable sites allowed to call a specimen MSI positive
msi_max_threshold=0.1


# Run the mSINGS analysis

BAMNAME=$(basename $BAM)
PFX=${BAMNAME%.*}

mkdir -p $OUTDIR/$PFX

echo
echo “Starting Analysis of $PFX” 
date +"%Y-%m-%d %H:%M" 

echo
echo "Making mpileup" 
date +"%Y-%m-%d %H:%M" 
samtools mpileup -f $REF_GENOME -d 100000 -A -E $BAM -l $INTERVALS_FILE | awk '{if($4 >= 6) print $0}' > $OUTDIR/$PFX/$PFX.mpileup  # make pileup if depth >= 6 reads

echo
echo "Varscan Readcounts start" 
date +"%Y-%m-%d %H:%M" 
java -Xmx4g -jar $VARSCAN readcounts $OUTDIR/$PFX/$PFX.mpileup --variants-file $INTERVALS_FILE --min-base-qual 10 --output-file $OUTDIR/$PFX/$PFX.msi_output &
wait

echo
echo "MSI Analyzer start"
date +"%Y-%m-%d %H:%M"

msi analyzer $OUTDIR/$PFX/$PFX.msi_output $BEDFILE -o $OUTDIR/$PFX/$PFX.msi.txt

echo
echo "MSI calls start" 
date +"%Y-%m-%d %H:%M" 
 
msi count_msi_samples $MSI_BASELINE $OUTDIR/$PFX -m $multiplier -t $msi_min_threshold $msi_max_threshold -o $OUTDIR/$PFX/$PFX.MSI_Analysis.txt

echo
echo “Completed Analysis of $PFX” 
date +"%Y-%m-%d %H:%M"  
echo








