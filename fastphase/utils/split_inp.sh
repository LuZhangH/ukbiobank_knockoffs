#!/bin/bash

##############################
# Parse input
##############################
PROGRAM_NAME=$0

function display_usage {
    echo "Usage: $PROGRAM_NAME -i input -a train -e test -n size"
    echo "  -i input           basename for INP/sample/legend files"
    echo "  -a train           basename for output train INP file"
    echo "  -e test            basename for output test INP file"
    echo "  -n size            size of test set"
    exit 1
}

# Default values for optional input arguments
RESET=0

# If less than 4 arguments supplied, display usage
if [  $# -le 5 ]
then
  display_usage
  exit 1
fi

# Parse arguments
echo "Parsed input arguments for "$PROGRAM_NAME":"
while getopts ":i:a:e:n:" opt; do
  case $opt in
    i)
      echo "  - input                 : $OPTARG" >&2
      INPUT=$OPTARG
      ;;
    a)
      echo "  - train                 : $OPTARG" >&2
      TRAIN=$OPTARG
      ;;
    e)
      echo "  - test                  : $OPTARG" >&2
      TEST=$OPTARG
      ;;
    n)
      echo "  - size of test set      : $OPTARG" >&2
      N_TEST=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done
echo ""

# Divide individuals
INDIVIDUALS=$INPUT".sample"

# Create list of lines
LINES_TRAIN=$TRAIN".lines"
LINES_TEST=$TEST".lines"
LINES_TMP=$TRAIN".lines.tmp"
INP_TRAIN=$TRAIN".inp"
INP_TEST=$TEST".inp"

# Divide individuals and create list of lines
SEED=2018
N_SAMPLES=$(awk 'END{print (NR-2)}' $INDIVIDUALS)
seq 1 $N_SAMPLES > $LINES_TMP
utils/shuffle.sh -s $SEED -n $N_TEST $LINES_TMP | sort -n > $LINES_TEST
comm -23 <( sort $LINES_TMP ) <( sort $LINES_TEST ) | sort -n > $LINES_TRAIN
rm $LINES_TMP

# Convert lists of lines into lists of lines for inp file
INP_LINES_TRAIN=$TRAIN".inp.lines"
INP_LINES_TEST=$TEST".inp.lines"
awk 'BEGIN{print 2; print 3}
    {print 4+3*($1-1); print 5+3*($1-1); print 6+3*($1-1)}' $LINES_TRAIN > $INP_LINES_TRAIN
awk 'BEGIN{print 2; print 3}
    {print 4+3*($1-1); print 5+3*($1-1); print 6+3*($1-1)}' $LINES_TEST > $INP_LINES_TEST

# Print first line of INP files
awk 'END{print NR}' $LINES_TRAIN > $INP_TRAIN   # Number of individuals in train set
awk 'END{print NR}' $LINES_TEST > $INP_TEST     # Number of individuals in test set

# Copy selected lines of INP files
awk -v linesfile=$INP_LINES_TRAIN -f utils/extract.awk $INPUT".inp" >> $INP_TRAIN
awk -v linesfile=$INP_LINES_TEST -f utils/extract.awk $INPUT".inp" >> $INP_TEST

# Convert lists of lines into lists of lines for sample file
SAMPLE_LINES_TRAIN=$TRAIN".sample.lines"
SAMPLE_LINES_TEST=$TEST".sample.lines"
awk 'BEGIN{print 1; print 2}
    {print 2+$1}' $LINES_TRAIN > $SAMPLE_LINES_TRAIN
awk 'BEGIN{print 1; print 2}
    {print 2+$1}' $LINES_TEST > $SAMPLE_LINES_TEST

# Copy selected lines of sample files
awk -v linesfile=$SAMPLE_LINES_TRAIN -f utils/extract.awk $INPUT".sample" > $TRAIN".sample"
awk -v linesfile=$SAMPLE_LINES_TEST -f utils/extract.awk $INPUT".sample" > $TEST".sample"

# Clean up temporary files
rm $LINES_TRAIN $LINES_TEST $INP_LINES_TRAIN $INP_LINES_TEST $SAMPLE_LINES_TRAIN $SAMPLE_LINES_TEST

# Copy legend files
cp $INPUT".legend" $TRAIN".legend"
cp $INPUT".legend" $TEST".legend"
