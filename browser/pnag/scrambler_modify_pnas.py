import sys
from Bio import SeqIO
from difflib import SequenceMatcher
import pandas as pd

seq_path = sys.argv[1]
res_path = sys.argv[2]
input_pna = sys.argv[3]

seqs = []
aso_seqs = []
output_df = pd.DataFrame(columns=["ASO", "ASO_seq", "SC_bases", "long_pur_stretch", "OT_TIR_0mm", "OT_TIR_1mm",
                                  "OT_TIR_2mm", "OT_TIR_3mm", "OT_tot_0mm", "OT_tot_1mm", "OT_tot_2mm", "OT_tot_3mm"])

# frst, check for input PNA sequence (fasta file) and add it to the output_df:
for record in SeqIO.parse(input_pna, "fasta"):
    aso_raw = record.seq

    aso_target = record.seq.reverse_complement()
    maxcomp_raw = SequenceMatcher(None, aso_raw, aso_target).find_longest_match(0, len(aso_raw), 0, len(aso_target)).size

    aso_name = "input PNA"
    # check for longest purine stretch and purine perc:
    pur = 0
    longest_purine_stretch_raw = 0
    curstretch = 0
    for base in aso_raw.__str__():
        if base in ["A", "G"]:
            curstretch += 1
            pur += 1
            if curstretch > longest_purine_stretch_raw:
                longest_purine_stretch_raw += 1
        else:
            curstretch = 0

    added_row = pd.Series([aso_name, aso_raw.__str__(), maxcomp_raw, longest_purine_stretch_raw, None, None, None, None,
                            None, None, None, None],
                          index=output_df.columns)
    # do this but using contat and not append: output_df = output_df.append(added_row, ignore_index=True)
    output_df = pd.concat([output_df, added_row.to_frame().transpose()], ignore_index=True)

    seqs.append(SeqIO.SeqRecord(aso_target, aso_name, description=""))
    aso_seqs.append(SeqIO.SeqRecord(aso_raw, aso_name, description=""))


for record in SeqIO.parse(seq_path, "fasta"):
    aso = record.seq

    aso_target = record.seq.reverse_complement()
    maxcomp = SequenceMatcher(None, aso, aso_target).find_longest_match(0, len(aso), 0, len(aso_target)).size
    # get identity between input PNA and current ASO:
    max_similarity_raw = SequenceMatcher(None, aso_raw, aso).find_longest_match(0, len(aso_raw), 0, len(aso)).size
    perc_similarity_raw = (max_similarity_raw / len(aso_raw)) * 100
    aso_name = record.id
    # check for longest purine stretch and purine perc:
    pur = 0
    longest_purine_stretch = 0
    curstretch = 0
    for base in aso.__str__():
        if base in ["A", "G"]:
            curstretch += 1
            pur += 1
            if curstretch > longest_purine_stretch:
                longest_purine_stretch += 1
        else:
            curstretch = 0

    # add to output_df. but ony if perc_similarity_raw < 50 and maxcomp_raw +- 1 of maxcomp and
    # longest_purine_stretch < longest_purine_stretch_raw+1
    if perc_similarity_raw < 30 and maxcomp_raw in [maxcomp, maxcomp-1, maxcomp+1] and \
            longest_purine_stretch < longest_purine_stretch_raw+1:
        added_row = pd.Series([aso_name, aso.__str__(), maxcomp, longest_purine_stretch, None, None, None, None, None,
                                 None, None, None],
                              index=output_df.columns)

        # do this but with concat! output_df = output_df.append(added_row, ignore_index=True):
        output_df = pd.concat([output_df, added_row.to_frame().transpose()], ignore_index=True)
        seqs.append(SeqIO.SeqRecord(aso_target, record.id, description=""))
        aso_seqs.append(SeqIO.SeqRecord(aso, record.id, description=""))
    else:
        print("ASO not added to output_df because of similarity to input PNA or purine stretch length.")
        print("perc_similarity_raw < 50: ", perc_similarity_raw)
        print("maxcomp_raw in [maxcomp, maxcomp-1, maxcomp+1]: ", maxcomp_raw, maxcomp)
        print("longest_purine_stretch < longest_purine_stretch_raw+1: ",longest_purine_stretch_raw, longest_purine_stretch)

SeqIO.write(seqs, res_path + "/reference_sequences/aso_targets.fasta", "fasta")
SeqIO.write(aso_seqs, res_path + "/reference_sequences/shuffled_sequences.fasta", "fasta")
output_df.to_csv(res_path + "/outputs/result_table.tsv", sep="\t", index=False)

print("finished modift_PNAs script!")

