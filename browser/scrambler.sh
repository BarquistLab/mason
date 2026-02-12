#!/bin/bash

# I start with assigning the flags (user inputs):


while getopts ":f:g:p:i:m:s:t:u:" flag; do
    case "${flag}" in
        f) fasta=${OPTARG} ;;
        g) gff=${OPTARG} ;;
        i) result_id=${OPTARG} ;;
        p) pna_input=${OPTARG} ;;
        m) mode=${OPTARG} ;;
        s) screen=${OPTARG} ;;
        t) target=${OPTARG} ;;
        u) use_ml=${OPTARG} ;;
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
  echo "start shuffling!" >> logfile_masonscript.log 2>&1
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

  if [ -n "$target" ]
  then
    echo "checker target mode: extracting target gene $target"

    # extract raw CDS transcripts (same as mason.sh)
    grep -P "\tCDS\t|\tsRNA\t|\tncRNA\t|\tgene\t" "$gff" |\
      awk -F'\t' 'BEGIN { OFS="\t" } {print $0}'| \
      awk -F'\t' 'BEGIN { OFS="\t" } { if ($9 ~ /locus_tag=/) { print $0 } else if ($9 ~ /ID=[^;]*/) { match($9, /ID=[^;]*/); print $0 ";locus_tag=" substr($9, RSTART+3, RLENGTH-3) } }' | \
      sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;locus_tag=([A-Za-z0-9_-]+).*)/\1\4\3/' | \
      sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;gene=([A-Za-z0-9_-]+).*)/\1\2;\4\3/' |
      grep ";locus_tag="> \
      "$REF/raw_transcripts_$GFF_NEW"

    # extract raw CDS FASTA for CAI calculation
    bedtools getfasta -s -fi "$fasta" -bed "$REF/raw_transcripts_$GFF_NEW" \
       -name+ -fo "$REF/raw_transcripts_$FASTA_NEW" >> logfile_masonscript.log 2>&1

    # calculate CAI
    cusp -sequence "$REF/raw_transcripts_$FASTA_NEW" -outfile "$REF/codon_occ.cusp" >> logfile_masonscript.log 2>&1
    grep -A 1 "$target" "$REF/raw_transcripts_$FASTA_NEW" | head -2 > "$REF/targetgene_fullseq.fasta"
    cai -seqall "$REF/targetgene_fullseq.fasta" -cfile "$REF/codon_occ.cusp" -outfile "$OUT/cai_value.txt" >> logfile_masonscript.log 2>&1
    sed -i 's/Sequence: [0-9]*-[0-9]*([+-]) CAI: //' "$OUT/cai_value.txt"

    # clean up raw transcript files
    rm -f "$REF/raw_transcripts_$FASTA_NEW" "$REF/raw_transcripts_$GFF_NEW" "$REF/codon_occ.cusp" "$REF/targetgene_fullseq.fasta"

    # extract target gene TIR regions from full transcripts
    grep -A 1 "$target" "$REF/full_transcripts_$FASTA_NEW" | \
      head -2 | sed -E 's/^([A-Z]{46}).*/\1/' > "$REF/targetgene_startreg.fasta"

    grep -A 1 "$target" "$REF/full_transcripts_$FASTA_NEW" | \
      head -2 | sed -E 's/^([A-Z]{60}).*/\1/' > "$REF/targetgene_mfe.fasta"

    # validate target gene was found in GFF
    if [ ! -s "$REF/targetgene_startreg.fasta" ] || [ "$(wc -l < "$REF/targetgene_startreg.fasta")" -lt 2 ]
    then
      echo "Target gene $target not found in GFF annotation!"
      echo "Target gene $target not found in GFF annotation!" >> "$RES/error.txt"
      touch "$RES/done.txt"
      exit 0
    fi

    # calculate MFE with RNAfold
    RNAfold -i "$REF/targetgene_mfe.fasta" --noPS > "$REF/intarna_output.fold"
    FOLD_FILE="$REF/intarna_output.fold"
    SEQ=$(awk 'NR==2' "$FOLD_FILE")
    STRUCT=$(awk 'NR==3 {print $1}' "$FOLD_FILE")
    MFE=$(grep -o '[(][[:space:]]*[-][0-9.]\+' "$FOLD_FILE" | tail -1 | tr -d '()[:space:]')
    echo "$MFE" > "$OUT/mfe_values.txt"

    # run checker_modify_pnas.py with target gene FASTA args
    python "./pnag/checker_modify_pnas.py" "$REF/PNA_sequence.fasta" "$RES" \
      "$REF/targetgene_startreg.fasta" "$REF/targetgene_mfe.fasta" >> logfile_masonscript.log 2>&1

    # validate that at least some ASOs target the gene
    if [ ! -s "$REF/aso_targets.fasta" ] || [ "$(wc -l < "$REF/aso_targets.fasta")" -lt 2 ]
    then
      echo "No ASOs target the specified gene $target!"
      echo "No ASOs target the specified gene $target! None of the provided ASO sequences match the translation initiation region." >> "$RES/error.txt"
      touch "$RES/done.txt"
      exit 0
    fi

    # Read SD position for VARNA highlighting (written by checker_modify_pnas.py)
    if [ -s "$OUT/sd_position.txt" ]; then
        read SD_START SD_END < "$OUT/sd_position.txt"
        SD_HIGHLIGHT="${SD_START}-${SD_END}:fill=#FFA500,outline=#FFA500,radius=10;"
    else
        SD_HIGHLIGHT=""
    fi
    VARNA_HIGHLIGHT="${SD_HIGHLIGHT}31-33:fill=#FF0000,outline=#FF0000,radius=10"

    # generate VARNA plots: one per-gene structure + per-ASO with binding highlighted
    # first: plain TIR structure (same as mason.sh)
    varna -sequenceDBN "$SEQ" \
      -structureDBN "$STRUCT" \
      -highlightRegion "$VARNA_HIGHLIGHT" \
      -title "Secondary structure of $target (MFE = $MFE kcal/mol)" \
      -titleSize 10 \
      -o "$OUT/varna_plot.svg"

    varna -sequenceDBN "$SEQ" \
      -structureDBN "$STRUCT" \
      -highlightRegion "$VARNA_HIGHLIGHT" \
      -title "Sec. structure of $target (MFE = $MFE kcal/mol)" \
      -titleSize 10 \
      -resolution "3.0" \
      -o "$OUT/varna_plot.png"

    # per-ASO VARNA plots: ASO binding drawn first with larger radius so it remains
    # visible behind the SD (orange) and start codon (red) regions
    if [ -f "$OUT/varna_positions.tsv" ]
    then
      tail -n +2 "$OUT/varna_positions.tsv" | while IFS=$'\t' read -r aso_name start_pos end_pos
      do
        if [ -n "$SD_HIGHLIGHT" ]; then
            HIGHLIGHT="$start_pos-$end_pos:fill=#4169E1,outline=#4169E1,radius=15;${SD_HIGHLIGHT}31-33:fill=#FF0000,outline=#FF0000,radius=10"
        else
            HIGHLIGHT="$start_pos-$end_pos:fill=#4169E1,outline=#4169E1,radius=15;31-33:fill=#FF0000,outline=#FF0000,radius=10"
        fi
        varna -sequenceDBN "$SEQ" \
          -structureDBN "$STRUCT" \
          -highlightRegion "$HIGHLIGHT" \
          -title "$aso_name binding on $target (MFE = $MFE kcal/mol)" \
          -titleSize 10 \
          -o "$OUT/varna_${aso_name}.svg"

        varna -sequenceDBN "$SEQ" \
          -structureDBN "$STRUCT" \
          -highlightRegion "$HIGHLIGHT" \
          -title "$aso_name on $target (MFE = $MFE kcal/mol)" \
          -titleSize 10 \
          -resolution "3.0" \
          -o "$OUT/varna_${aso_name}.png"
      done
    fi
  else
    python "./pnag/checker_modify_pnas.py" "$REF/PNA_sequence.fasta" "$RES" >> logfile_masonscript.log 2>&1
  fi
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
  if [ -n "$target" ]
  then
    Rscript ./pnag/summarize_ots_checker.R "$OUT" "$screen" "$target" "$use_ml" >> logfile_masonscript.log 2>&1
    if [ "$use_ml" = "yes" ]; then
      python ./pnag/ML_run.py "$OUT" >> logfile_masonscript.log 2>&1
    fi
    Rscript ./pnag/make_final_table_mason.R "$OUT" "$screen" "$use_ml" >> logfile_masonscript.log 2>&1
  else
    Rscript ./pnag/summarize_ots_checker.R "$OUT" "$screen"  >> logfile_masonscript.log 2>&1
  fi
else
  echo "Somethings wrong!"
fi

touch "$RES/$target"

if [ -n "$target" ] && [ "$mode" == "checker" ]
then
  cleanup_intermediate_files "checker_target"
else
  cleanup_intermediate_files "scrambler"
fi

echo "MASON finished" >> logfile_masonscript.log
