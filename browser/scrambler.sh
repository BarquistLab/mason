#!/bin/bash

# I start with assigning the flags (user inputs):


while getopts ":f:g:p:i:m:s:" flag; do
    case "${flag}" in
        f) fasta=${OPTARG} ;;
        g) gff=${OPTARG} ;;
        i) result_id=${OPTARG} ;;
        p) pna_input=${OPTARG} ;;
        m) mode=${OPTARG} ;;
        s) screen=${OPTARG} ;;
        \?) echo "Invalid option: -$OPTARG" >&2
            exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2
            exit 1 ;;
    esac
done

# re-activate conda environment:
. /home/jakob/miniconda3/etc/profile.d/conda.sh
conda activate browser

# Source shared pipeline functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common_pipeline.sh"


# Setup directories using shared function
setup_directories

# I print them out to be sure it worked out:
echo "fasta: $fasta";
echo "gff: $gff";
# change gff and fasta names to replace spaces with underscores (for mason, because mason doesn't like spaces in file names)
GFF_NEW="$(echo "$GFF" | tr ' ' '_')"
FASTA_NEW="$(echo "$FASTA" | tr ' ' '_')"

echo "$GFF $GFF_NEW"
echo "$FASTA $FASTA_NEW"

# Copy reference files using shared function (scrambler doesn't rename, so GFF_NEW/FASTA_NEW not set)
copy_reference_files

# extract the full transcripts from the gff (-30/+30):
extract_gff_transcripts "$gff" "$REF/full_transcripts_$GFF_NEW"

# extract the gene lengths:
extract_gene_lengths

# change the gff entries that go too far:
Rscript pnag/modify_gff.R "$REF/full_transcripts_$GFF_NEW" "$REF/genelengths.tsv" "$WARNINGS"
echo "adjusted gff successfully"

echo "start running bedtools"
# I extract the fasta files from the gff using bedtools:
bedtools getfasta -s -fi "$fasta" -bed "$REF/full_transcripts_$GFF_NEW"  \
	 -name+ -fo "$REF/full_transcripts_$FASTA_NEW"

echo "PNA $pna_input put in"


# shuffle if there is no input "checker" in the -m flag

if [ -z "$mode" ]
then
  echo ">PNA" > "$REF/PNA_sequence.fasta"
  echo "$pna_input" >> "$REF/PNA_sequence.fasta"
  echo "start shuffling!" >> logfile_manuscript.log >> logfile_masonscript.log 2>&1
  esl-shuffle -N 500 -o "$REF/shuffled_sequences.fasta" "$REF/PNA_sequence.fasta"
  # use sed to change all -shuffled- to _scr_
  sed -i 's/-shuffled-/_scr_/g' "$REF/shuffled_sequences.fasta"
  # add a 0 to all single digit numbers:
  sed -Ei 's/_scr_([0-9])$/_scr_0\1/g' "$REF/shuffled_sequences.fasta"
  # add a 0 to all double digit numbers:
  sed -Ei 's/_scr_([0-9][0-9])$/_scr_0\1/g' "$REF/shuffled_sequences.fasta"

  echo "modify PNAS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  python "./pnag/scrambler_modify_pnas.py" "$REF/shuffled_sequences.fasta" "$RES" "$REF/PNA_sequence.fasta" >> logfile_masonscript.log 2>&1

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
  # if input $pna_input is a string starting with >
  if [[ $pna_input == ">"* ]]
  then
      echo "$pna_input" > "$REF/PNA_sequence.fasta"
    echo "modify PNAS!!"

  else
    echo ">PNA" > "$REF/PNA_sequence.fasta"
    echo "$pna_input" >> "$REF/PNA_sequence.fasta"
    # else write it to PNA_sequence.fasta with >PNA as header
  fi
  python "./pnag/checker_modify_pnas.py" "$REF/PNA_sequence.fasta" "$RES" #>> logfile_masonscript.log 2>&1
else
  echo "Somethings wrong!"
fi



run_seqmap_full_transcriptome
run_optional_screening

# Process mismatches using shared function (offset=-32 for scrambler)
run_seqmap_and_process_mismatches -32

echo "summarize off-targets"

if [ -z "$mode" ]
then
  python ./pnag/summarize_ots_scrambler.py "$OUT" "$screen"  >> logfile_masonscript.log 2>&1
  Rscript ./pnag/plot_ots_scrambler.R "$OUT" "$screen" >> logfile_masonscript.log 2>&1
elif [ "$mode" == "checker" ]
then
  Rscript ./pnag/summarize_ots_checker.R "$OUT" "$screen"  >> logfile_masonscript.log 2>&1
else
  echo "Somethings wrong!"
fi

touch "$RES/$target"

echo "MASON finished" >> logfile_masonscript.log
