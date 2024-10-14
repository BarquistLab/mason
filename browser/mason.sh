#!/bin/bash

# I start with assigning the flags (user inputs):


while getopts f:g:t:l:m:p:i:b:s: flag
do
    # shellcheck disable=SC2220
    case "${flag}" in
	f) fasta=${OPTARG};;
	g) gff=${OPTARG};;
	t) target=${OPTARG};;
	l) length=${OPTARG};;
	m) mismatches=${OPTARG};;
	i) result_id=${OPTARG};;
	p) pna_input=${OPTARG};;
	b) bases_before=${OPTARG};;
	s) screen=${OPTARG};;
    esac
done

# re-activate conda environment:
. /home/jakob/miniconda3/etc/profile.d/conda.sh
conda activate browser


# I print them out to be sure it worked out:
echo "bases_before= $bases_before"
echo "fasta: $fasta";
echo "gff: $gff";
echo "screen: $screen";

mkdir "./pnag/static/data/$result_id"

RES="./pnag/static/data/$result_id"
REF="$RES/reference_sequences"
PRESETS="./pnag/static/data/presets"
OUT="$RES/outputs"
GFF="$(basename -- "$gff")"
FASTA="$(basename -- "$fasta")"  # get base names of files wo paths
WARNINGS="$RES/warnings.txt"


echo "$REF"
echo "$OUT"
GFF_NEW="$(echo "$GFF" | tr ' ' '_')"
FASTA_NEW="$(echo "$FASTA" | tr ' ' '_')"

echo "$GFF $GFF_NEW"
echo "$FASTA $FASTA_NEW"

find ./pnag/static/data/""20*  -mtime +30  -delete

mkdir -p $REF $OUT
touch "$WARNINGS"

scp "$gff" "$REF/$GFF_NEW"
scp "$fasta" "$REF/$FASTA_NEW"

gff="$REF/$GFF_NEW"
fasta="$REF/$FASTA_NEW"



# same for full regions (change to whole CDS and 30 nt upstream):
grep -P "\tCDS\t|\tsRNA\t|\tncRNA\t|\tgene\t" $gff |\
        awk -F'\t' 'BEGIN { OFS="\t" } {if ($7=="-") {$5=$5+30} else { $4=$4-30} print $0}'| \
    sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;locus_tag=([A-Za-z0-9_]+).*)/\1\4\3/' | \
    sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;gene=([A-Za-z0-9_]+).*)/\1\2;\4\3/' |
    grep ";locus_tag="> \
    "$REF/full_transcripts_$GFF_NEW"

bioawk -c fastx '{ print $name, length($seq) }' < "$fasta"  > "$REF/genelengths.tsv"

# change the gff entries that go too far:
Rscript pnag/modify_gff.R "$REF/full_transcripts_$GFF_NEW" "$REF/genelengths.tsv" "$WARNINGS"



echo "start running bedtools"
# I extract the fasta files from the gff using bedtools:
bedtools getfasta -s -fi "$fasta" -bed "$REF/full_transcripts_$GFF_NEW"  \
	 -name+ -fo "$REF/full_transcripts_$FASTA_NEW" >> logfile_masonscript.log 2>&1

echo "starting if statement"
echo "$pna_input"


if [ -z  "$pna_input" ];
then
    echo "no PNA put in"
    # Now I create a list of all PNAs:
    grep -A 1 $target "$REF/full_transcripts_$FASTA_NEW" | \
	sed -E 's/^([A-Z]{46}).*/\1/' > "$REF/targetgene_startreg.fasta"   |# select -30 to + 16 region
    # Now I run the python script which I wrote to design PNAs:
    echo "$length"
    echo "$result_id"
    python ./pnag/make_pnas.py $length $RES $bases_before >> logfile_masonscript.log 2>&1

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


#Now I run seqmap on start regions and whole transcriptome:
echo "run seqmap"
seqmap "$mismatches" "$REF/aso_targets.fasta" "$REF/full_transcripts_$FASTA_NEW" \
       "$OUT/offtargets_fulltranscripts.tab" /output_all_matches \
       /forward_strand /output_statistics /available_memory:5000 >> logfile_masonscript.log 2>&1



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
	     #$2=$2-32
	     print $0
	}' |  sed -E 's/^([^;:]*)::/\1;\1::/'| \
	    sed -E 's/^([^;:]*);([^:]*)::[^\(]*\(([\+\-])\)/\1\t\2\t\3/'| \
	     sed -E 's/^([A-Z][A-Z]_[^ ]+) ([^\t]+)/\1\t\2\tU/'| \
	     sed -E 's/^([A-Z][A-Z]_[^_]+)_([^\(]+)\(([\+\-])\)/\1\t\2\t\3/'  >> "$NEWNAME"
    rm "$NAME"
    
done

rm -rf $OUT/*_sorted_sorted.tab
python ./pnag/summarize_offtargets.py "$OUT" "$screen" >> logfile_masonscript.log 2>&1

touch "$RES/$target"

rm -rf "$OUT/offtargets_fulltranscripts_sorted.tab"
rm -rf "$OUT/offtargets_startregions_sorted.csv"

echo "MASON finished" >> logfile_masonscript.log 

