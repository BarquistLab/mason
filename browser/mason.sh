#!/bin/bash

# I start with assigning the flags (user inputs):


while getopts f:g:t:l:m:p:i:b:s:u: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
	f) fasta=${OPTARG};;
	g) gff=${OPTARG};;
	t) target=${OPTARG};;
	l) length=${OPTARG};;
	i) result_id=${OPTARG};;
	p) pna_input=${OPTARG};;
	b) bases_before=${OPTARG};;
	s) screen=${OPTARG};;
	u) use_ml=${OPTARG};;
    esac
done

# re-activate conda environment:
. /home/jakob/miniconda3/etc/profile.d/conda.sh
conda activate browser
export HOME=${HOME:-/home/jakob}

# Source shared pipeline functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common_pipeline.sh"

# if bases_before is not set, set it to length - 3
if [ -z "$bases_before" ]; then
    bases_before=$((0))
fi

# I print them out to be sure it worked out:
echo "bases_before= $bases_before"
echo "fasta: $fasta";
echo "gff: $gff";
echo "screen: $screen";

# Setup directories using shared function
setup_directories

# Copy essential genes list if present (one level up from gene-specific RES)
if [[ -f "$RES/../essential_genes.txt" ]]; then
    cp "$RES/../essential_genes.txt" "$OUT/essential_genes.txt"
fi

echo "$REF"
echo "$OUT"
# change gff and fasta names to replace spaces with underscores (for mason, because mason doesn't like spaces in file names)
GFF_NEW="$(echo "$GFF" | tr ' ' '_')"
FASTA_NEW="$(echo "$FASTA" | tr ' ' '_')"

echo "$GFF $GFF_NEW"
echo "$FASTA $FASTA_NEW"

# Copy reference files using shared function (uses GFF_NEW/FASTA_NEW for mason)
copy_reference_files


# extract full regions (change to whole CDS and 30 nt upstream):
extract_gff_transcripts "$gff" "$REF/full_transcripts_$GFF_NEW"

# extract only raw cds/genes, without -30/+30:
grep -P "\tCDS\t|\tsRNA\t|\tncRNA\t|\tgene\t" $gff |\
      awk -F'\t' 'BEGIN { OFS="\t" } {print $0}'| \
        # if theres no locus_tag, extract ID=(...) and add the extracted ID as locus_tag=... in the end of the line. use awk
        awk -F'\t' 'BEGIN { OFS="\t" } { if ($9 ~ /locus_tag=/) { print $0 } else if ($9 ~ /ID=[^;]*/) { match($9, /ID=[^;]*/); print $0 ";locus_tag=" substr($9, RSTART+3, RLENGTH-3) } }' | \
        sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;locus_tag=([A-Za-z0-9_-]+).*)/\1\4\3/' | \
        sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;gene=([A-Za-z0-9_-]+).*)/\1\2;\4\3/' |
    grep ";locus_tag="> \
    "$REF/raw_transcripts_$GFF_NEW"

extract_gene_lengths

# change the gff entries that go too far:
Rscript pnag/modify_gff.R "$REF/full_transcripts_$GFF_NEW" "$REF/genelengths.tsv" "$WARNINGS" "$REF/raw_transcripts_$GFF_NEW"



echo "start running bedtools"
# I extract the fasta files from the gff using bedtools:
bedtools getfasta -s -fi "$fasta" -bed "$REF/full_transcripts_$GFF_NEW"  \
	 -name+ -fo "$REF/full_transcripts_$FASTA_NEW" >> logfile_masonscript.log 2>&1

echo "start calculation of CAI"

# same for raw transcripts:
bedtools getfasta -s -fi "$fasta" -bed "$REF/raw_transcripts_$GFF_NEW"  \
   -name+ -fo "$REF/raw_transcripts_$FASTA_NEW" >> logfile_masonscript.log 2>&1

# I calculate the CAI values for the raw transcripts:
cusp -sequence "$REF/raw_transcripts_$FASTA_NEW" -outfile "$REF/codon_occ.cusp" >> logfile_masonscript.log 2>&1
# calculate CAI for target gene:
# first, create gene fasta
grep -A 1 "$target" "$REF/raw_transcripts_$FASTA_NEW" | \
  # make sure only first 2 lines are printed
  head -2 > "$REF/targetgene_fullseq.fasta"
# calculate CAI
cai -seqall "$REF/targetgene_fullseq.fasta" -cfile "$REF/codon_occ.cusp" -outfile "$OUT/cai_value.txt" >> logfile_masonscript.log 2>&1
# only preserve cai (now its: Sequence: 205125-208608(+) CAI: 0.726)
sed -i 's/Sequence: [0-9]*-[0-9]*([+-]) CAI: //' "$OUT/cai_value.txt"

# delete all files that are not needed anymore
rm "$REF/raw_transcripts_$FASTA_NEW"

rm "$REF/raw_transcripts_$GFF_NEW"
rm "$REF/codon_occ.cusp"
rm "$REF/targetgene_fullseq.fasta"
echo "finished CAI calculation"


echo "starting if statement"
echo "$pna_input"


if [ -z  "$pna_input" ];
then
    echo "no PNA put in"
    # Now I create a list of all PNAs:
    grep -A 1 $target "$REF/full_transcripts_$FASTA_NEW" | \
	    sed -E 's/^([A-Z]{46}).*/\1/' > "$REF/targetgene_startreg.fasta" # select -30 to + 16 region
	  # do same but get -30 to +30 region
	  grep -A 1 $target "$REF/full_transcripts_$FASTA_NEW" | \
      sed -E 's/^([A-Z]{60}).*/\1/' > "$REF/targetgene_mfe.fasta"

    # calculate MFE for target gene start region w rnafold
    RNAfold -i "$REF/targetgene_mfe.fasta" --noPS > "$REF/intarna_output.fold" #>> logfile_masonscript.log 2>&1

    FOLD_FILE="$REF/intarna_output.fold"
    SEQ=$(awk 'NR==2' "$FOLD_FILE")
    STRUCT=$(awk 'NR==3 {print $1}' "$FOLD_FILE")
    MFE=$(grep -o '[(][[:space:]]*[-][0-9.]\+' "$FOLD_FILE" | tail -1 | tr -d '()[:space:]')
    # save MFE value as mfe_values.txt
    echo "$MFE" > "$OUT/mfe_values.txt"

    # Now I run the python script which I wrote to design PNAs:
    echo "$length"
    echo "$result_id"
    python ./pnag/make_pnas.py "$length" "$RES" "$bases_before" >> logfile_masonscript.log 2>&1

    # Read SD position for VARNA highlighting (written by make_pnas.py)
    if [ -s "$OUT/sd_position.txt" ]; then
        read SD_START SD_END < "$OUT/sd_position.txt"
        SD_HIGHLIGHT="${SD_START}-${SD_END}:fill=#FFA500,outline=#FFA500,radius=10;"
    else
        SD_HIGHLIGHT=""
    fi
    VARNA_HIGHLIGHT="${SD_HIGHLIGHT}31-33:fill=#FF0000,outline=#FF0000,radius=10"

    echo "Running VARNA SVG..." >> logfile_masonscript.log
    xvfb-run -a /home/jakob/bin/varna -sequenceDBN "$SEQ" \
      -structureDBN "$STRUCT" \
      -highlightRegion "$VARNA_HIGHLIGHT" \
      -title "Secondary structure of $target (MFE = $MFE kcal/mol)" \
      -titleSize 10 \
      -o "$OUT/varna_plot.svg" >> logfile_masonscript.log 2>&1

    # same but output a png file
    echo "Running VARNA PNG..." >> logfile_masonscript.log
    xvfb-run -a /home/jakob/bin/varna -sequenceDBN "$SEQ" \
      -structureDBN "$STRUCT" \
      -highlightRegion "$VARNA_HIGHLIGHT" \
      -title "Sec. structure of $target (MFE = $MFE kcal/mol)" \
      -titleSize 10 \
      -resolution "3.0" \
      -o "$OUT/varna_plot.png" >> logfile_masonscript.log 2>&1

else
    echo "PNA $pna_input put in"
    python ./pnag/modify_PNAs.py "$pna_input" >> logfile_masonscript.log 2>&1
fi

if [ ! -s "$REF/aso_targets.fasta" ]
then
  echo "No PNAs created!"
  touch "$RES/../done.txt"
  echo  "$target" >> "$RES/../error.txt"
fi


run_seqmap_full_transcriptome
run_optional_screening

# Process mismatches using shared function (offset=0 for mason)
run_seqmap_and_process_mismatches 0

Rscript ./pnag/summarize_ots.R "$OUT" "$target" "$screen" >> logfile_masonscript.log 2>&1

if [ "$use_ml" = "yes" ]; then
    python ./pnag/ML_run.py "$OUT" >> logfile_masonscript.log 2>&1
fi

Rscript ./pnag/make_final_table_mason.R "$OUT" "$screen" "$use_ml" >> logfile_masonscript.log 2>&1

touch "$RES/$target"

cleanup_intermediate_files "mason"

echo "MASON finished" >> logfile_masonscript.log
