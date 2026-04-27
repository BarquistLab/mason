#!/bin/bash
# common_pipeline.sh — shared functions for mason.sh and scrambler.sh

setup_directories() {
    # Creates RES, REF, OUT directories and sets common variables.
    # Expects: $result_id, $gff, $fasta to be set
    # Sets: RES, REF, PRESETS, OUT, GFF, FASTA, WARNINGS

    mkdir -p "./pnag/static/data/$result_id"

    RES="./pnag/static/data/$result_id"
    REF="$RES/reference_sequences"
    PRESETS="./pnag/static/data/presets"
    OUT="$RES/outputs"
    GFF="$(basename -- "$gff")"
    FASTA="$(basename -- "$fasta")"
    WARNINGS="$RES/warnings.txt"

    # delete old files (>30 days)
    find ./pnag/static/data/20*  -mtime +30  -delete 2>/dev/null

    mkdir -p "$REF" "$OUT"
    touch "$WARNINGS"
}

copy_reference_files() {
    # Copy FASTA/GFF to reference directory.
    # Expects: $gff, $fasta, $REF, $GFF, $FASTA to be set
    # Optional: $GFF_NEW, $FASTA_NEW for renamed copies (mason uses space-replaced names)
    local gff_dest="${GFF_NEW:-$GFF}"
    local fasta_dest="${FASTA_NEW:-$FASTA}"

    scp "$gff" "$REF/$gff_dest"
    scp "$fasta" "$REF/$fasta_dest"

    # Update gff/fasta to point to copied files
    gff="$REF/$gff_dest"
    fasta="$REF/$fasta_dest"
}

extract_gff_transcripts() {
    # The 6-line grep/awk/sed pipeline to extract full transcripts from GFF.
    # Expects: $gff, $REF, $GFF_BASENAME (the basename to use for output files)
    local gff_src="$1"
    local output_file="$2"

    grep -P "\tCDS\t|\tsRNA\t|\tncRNA\t|\tgene\t" "$gff_src" |\
        awk -F'\t' 'BEGIN { OFS="\t" } {if ($7=="-") {$5=$5+30} else { $4=$4-30} print $0}'| \
        awk -F'\t' 'BEGIN { OFS="\t" } { if ($9 ~ /locus_tag=/) { print $0 } else if ($9 ~ /ID=[^;]*/) { match($9, /ID=[^;]*/); print $0 ";locus_tag=" substr($9, RSTART+3, RLENGTH-3) } }' | \
        sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;locus_tag=([A-Za-z0-9_-]+).*)/\1\4\3/' | \
        sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;gene=([A-Za-z0-9_-]+).*)/\1\2;\4\3/' |
        grep ";locus_tag=" > \
        "$output_file"
}

run_seqmap_and_process_mismatches() {
    # Run seqmap on off-targets and process mismatch positions with AWK.
    # Args:
    #   $1 — coordinate_offset (0 for mason, -32 for scrambler)
    #   Remaining seqmap invocations should be done before calling this function.
    # Expects: $OUT to be set, offtargets*.tab files to exist in $OUT
    local coord_offset="$1"

    echo "determine mismatch positions"
    for NAME in "$OUT"/offtargets*.tab
    do
        echo "$NAME"
        NEWNAME=${NAME%.tab}_sorted.tab
        head -1 "$NAME" | sed -E "s/(.*)/\\1\tmismatch_positions\tlongest_stretch\tbinding_sequence/" |
        sed -E  's/^trans_id/locus_tag\tgene_name\tstrand/' > "$NEWNAME"

        echo "$NEWNAME"
        sed 1d "$NAME" |\
        sort -t "$(printf "\t")"  -k4 |\
        awk -F'\t' -v offset="$coord_offset" 'BEGIN {OFS="\t"; pos=0}
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
             if (offset != 0) { $2=$2+offset }
             print $0
        }' |  sed -E 's/^([^;:]*)::/\1;\1::/'| \
            sed -E 's/^([^;:]*);([^:]*)::[^\(]*\(([\+\-])\)/\1\t\2\t\3/'| \
             sed -E 's/^([A-Z][A-Z]_[^ ]+) ([^\t]+)/\1\t\2\tU/'| \
             sed -E 's/^([A-Z][A-Z]_[^_]+)_([^\(]+)\(([\+\-])\)/\1\t\2\t\3/'  >> "$NEWNAME"
        rm "$NAME"

    done

    rm -rf "$OUT"/*_sorted_sorted.tab
}

extract_gene_lengths() {
    # Extract gene lengths from FASTA using bioawk.
    # Expects: $fasta, $REF to be set
    bioawk -c fastx '{ print $name, length($seq) }' < "$fasta" > "$REF/genelengths.tsv"
}

run_seqmap_full_transcriptome() {
    # Run seqmap with 3 mismatches on full transcriptome.
    # Expects: $REF, $FASTA_NEW, $OUT to be set
    echo "run seqmap"
    seqmap 3 "$REF/aso_targets.fasta" "$REF/full_transcripts_$FASTA_NEW" \
           "$OUT/offtargets_fulltranscripts.tab" /output_all_matches \
           /forward_strand /output_statistics /available_memory:5000 >> logfile_masonscript.log 2>&1
}

run_optional_screening() {
    # Run seqmap for HMP microbiome or human transcriptome screening (0 mismatches).
    # Expects: $screen, $REF, $PRESETS, $OUT to be set
    if [[ $screen = "microbiome" ]]
    then
        seqmap 0 "$REF/aso_targets.fasta" "$PRESETS/start_regions_HMP.fasta" \
           "$OUT/offtargets_microbiome.tab" /output_all_matches \
           /forward_strand /output_statistics /available_memory:5000 >> logfile_masonscript.log 2>&1
    elif [[ $screen = "human" ]]
    then
        seqmap 0 "$REF/aso_targets.fasta" "$PRESETS/GRCh38_latest_rna.fna" \
           "$OUT/offtargets_human.tab" /output_all_matches \
           /forward_strand /output_statistics /available_memory:5000 >> logfile_masonscript.log 2>&1
    fi
}

cleanup_intermediate_files() {
    # Remove intermediate files that aren't downloadable from the web UI.
    # Args: $1 — mode ("mason", "scrambler", or "checker_target")
    # Expects: $REF, $OUT to be set
    local mode="$1"

    # Common reference_sequences/ files
    rm -f "$REF"/full_transcripts_*
    rm -f "$REF/genelengths.tsv"
    rm -f "$REF/aso_targets.fasta"
    rm -f "$REF"/*.fai

    # Mode-specific reference_sequences/ files
    if [[ "$mode" == "mason" || "$mode" == "checker_target" ]]; then
        rm -f "$REF/targetgene_startreg.fasta"
        rm -f "$REF/targetgene_mfe.fasta"
        rm -f "$REF/intarna_output.fold"
        rm -f "$REF/dot.ps"
    fi
    if [[ "$mode" == "scrambler" || "$mode" == "checker_target" ]]; then
        rm -f "$REF/PNA_sequence.fasta"
        rm -f "$REF/dot.ps"
    fi

    # Common outputs/ files
    rm -f "$OUT/df_plot.csv"
    rm -f "$OUT/result_table.tsv"
    rm -f "$OUT"/*_sorted.tab
    rm -f "$OUT/essential_genes.txt"

    # Mason and checker_target outputs/ files
    if [[ "$mode" == "mason" || "$mode" == "checker_target" ]]; then
        rm -f "$OUT/saved_table_ml.csv"
        rm -f "$OUT/cai_value.txt"
        rm -f "$OUT/mfe_values.txt"
    fi
}
