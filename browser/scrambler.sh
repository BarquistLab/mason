#!/bin/bash

# I start with assigning the flags (user inputs):


while getopts ":f:g:p:i:m:" flag; do
    case "${flag}" in
        f) fasta=${OPTARG} ;;
        g) gff=${OPTARG} ;;
        i) result_id=${OPTARG} ;;
        p) pna_input=${OPTARG} ;;
        m) mode=${OPTARG} ;;
        \?) echo "Invalid option: -$OPTARG" >&2
            exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2
            exit 1 ;;
    esac
done

# re-activate conda environment:
. /home/jakob/miniconda3/etc/profile.d/conda.sh
conda activate browser


# I print them out to be sure it worked out:
echo "fasta: $fasta";
echo "gff: $gff";

mkdir "./pnag/static/data/$result_id"

RES="./pnag/static/data/$result_id"
REF="$RES/reference_sequences"
PRESETS="./pnag/static/data/presets"
OUT="$RES/outputs"
GFF="$(basename -- $gff)"
FASTA="$(basename -- $fasta)"  # get base names of files wo paths
WARNINGS="$RES/warnings.txt"

# delete old files:
find ./pnag/static/data/20*  -mtime +30  -delete

# make directories & create warnings file:
mkdir -p $REF $OUT
touch "$WARNINGS"

# copy files to reference directory:
scp "$gff" "$REF/$GFF"
scp "$fasta" "$REF/$FASTA"

# extract the full transcripts from the gff (-30/+30):
grep -P "\tCDS\t|\tsRNA\t|\tncRNA\t|\tgene\t" $gff |\
        awk -F'\t' 'BEGIN { OFS="\t" } {if ($7=="-") {$5=$5+30} else { $4=$4-30} print $0}'| \
    sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;locus_tag=([^;]+).*)/\1\4\3/' | \
    sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;gene=([^;]+).*)/\1\2;\4\3/' |
    grep ";locus_tag="> \
    "$REF/full_transcripts_$GFF"

# extract the gene lengths:
bioawk -c fastx '{ print $name, length($seq) }' < "$fasta"  > "$REF/genelengths.tsv"

# change the gff entries that go too far:
Rscript pnag/modify_gff.R "$REF/full_transcripts_$GFF" "$REF/genelengths.tsv" "$WARNINGS"
echo "adjusted gff successfully"

echo "start running bedtools"
# I extract the fasta files from the gff using bedtools:
bedtools getfasta -s -fi $fasta -bed "$REF/full_transcripts_$GFF"  \
	 -name+ -fo "$REF/full_transcripts_$FASTA"

echo "PNA $pna_input put in"


# shuffle if there is no input "checker" in the -m flag

if [ -z "$mode" ]
then
  echo ">PNA" > "$REF/PNA_sequence.fasta"
  echo "$pna_input" >> "$REF/PNA_sequence.fasta"
  echo "start shuffling!"
  esl-shuffle -N 500 -o "$REF/shuffled_sequences.fasta" "$REF/PNA_sequence.fasta"
  # use sed to change all -shuffled- to _scr_
  sed -i 's/-shuffled-/_scr_/g' "$REF/shuffled_sequences.fasta"
  # add a 0 to all single digit numbers:
  sed -Ei 's/_scr_([0-9])$/_scr_0\1/g' "$REF/shuffled_sequences.fasta"
  # add a 0 to all double digit numbers:
  sed -Ei 's/_scr_([0-9][0-9])$/_scr_0\1/g' "$REF/shuffled_sequences.fasta"

  echo "modify PNAS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  python "./pnag/scrambler_modify_pnas.py" "$REF/shuffled_sequences.fasta" "$RES" "$REF/PNA_sequence.fasta" #>> logfile_masonscript.log 2>&1

  # check if aso_targets.fasta has less than 3 lines:
  if [ "$(wc -l < "$REF/aso_targets.fasta")" -lt 3 ]
  then
    echo "No PNAs created!"
    touch "$RES/../done.txt"
    echo  "No PNAs created!" >> "$RES/../error.txt"
  fi
elif [ "$mode" == "checker" ]
then
  echo "no shuffling!"
  echo "$pna_input" > "$REF/PNA_sequence.fasta"
  echo "modify PNAS!!"
  python "./pnag/checker_modify_pnas.py" "$REF/PNA_sequence.fasta" "$RES" #>> logfile_masonscript.log 2>&1
else
  echo "Somethings wrong!"
fi





#Now I run seqmap on start regions and whole transcriptome:
echo "run seqmap"
seqmap 3 "$REF/aso_targets.fasta" "$REF/full_transcripts_$FASTA" \
       "$OUT/offtargets_fulltranscripts.tab" /output_all_matches \
       /forward_strand /output_statistics /available_memory:5000 >> logfile_masonscript.log 2>&1


echo "determine mismatch positions"
# I use awk to determine the mismatch positions:
for NAME in "$OUT"/offtargets*.tab
do
    echo "$NAME"
    NEWNAME=${NAME%.tab}_sorted.tab
    head -1 $NAME | sed -E "s/(.*)/\\1\tmismatch_positions\tlongest_stretch\tbinding_sequence/" |
    sed -E  's/^trans_id/locus_tag\tgene_name\tstrand/' > $NEWNAME

    echo "$NEWNAME"
    sed 1d $NAME |\
	sort -t "$(printf "\t")"  -k4 |\
	awk -F'\t' 'BEGIN {OFS="\t"; pos=0}
	{
	    max=length($3)
	    mm="none"
	    stretch=0
	    longest_stretch=0
	    binding_sequence=""
	    longest_binding_sequence=""
	     for(i=1; pos == 0 && i <= max; i++)
	     {
		    v1=substr($3, i, 1)
		    v2=substr($5, i, 1)
		    if(v1 != v2)
		     {
		       stretch=0
		       binding_sequence=""
		       if(mm=="none") {mm=i} else {mm=mm ";" i}
		     }
		     else
		     {
		      binding_sequence=binding_sequence v1
          stretch++
          if(stretch > longest_stretch)
          {
          longest_stretch=stretch
          longest_binding_sequence=binding_sequence
          }
		     }
	     }
	     $7=mm
	     $8=longest_stretch
	     $9=longest_binding_sequence
	     $2=$2-32
	     print $0
	}' |  sed -E 's/^([^;:]*)::/\1;\1::/'| \
	    sed -E 's/^([^;:]*);([^:]*)::[^\(]*\(([\+\-])\)/\1\t\2\t\3/'| \
	     sed -E 's/^([A-Z][A-Z]_[^ ]+) ([^\t]+)/\1\t\2\tU/'| \
	     sed -E 's/^([A-Z][A-Z]_[^_]+)_([^\(]+)\(([\+\-])\)/\1\t\2\t\3/'  >> "$NEWNAME"
    rm "$NAME"

done

rm -rf $OUT/*_sorted_sorted.tab

echo "summarize off-targets"

if [ -z "$mode" ]
then
  python ./pnag/summarize_ots_scrambler.py "$OUT" # >> logfile_masonscript.log 2>&1
elif [ "$mode" == "checker" ]
then
  Rscript ./pnag/summarize_ots_checker.R "$OUT" # >> logfile_masonscript.log 2>&1
else
  echo "Somethings wrong!"
fi

touch "$RES/$target"
# remove offtargets_fulltranscripts_sorted.tab
rm -rf "$OUT/offtargets_fulltranscripts_sorted.tab"
rm -rf "$OUT/offtargets_startregions_sorted.csv"

echo "MASON finished" >> logfile_masonscript.log

