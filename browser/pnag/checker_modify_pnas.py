import sys
from Bio import SeqIO
import pandas as pd
from pna_utils import calculate_purine_stats, calculate_self_complementarity

res_path = sys.argv[2]
input_pna = sys.argv[1]

seqs = []
aso_seqs = []
output_df = pd.DataFrame(columns=["ASO", "ASO_seq", "SC_bases", "pur_perc", "long_pur_stretch", "OT_TIR_0mm", "OT_TIR_1mm",
                                  "OT_TIR_2mm", "OT_TIR_3mm", "OT_tot_0mm", "OT_tot_1mm", "OT_tot_2mm", "OT_tot_3mm"])

# frst, check for input PNA sequence (fasta file) and add it to the output_df:
for record in SeqIO.parse(input_pna, "fasta"):
    aso_raw = record.seq

    maxcomp_raw, aso_target = calculate_self_complementarity(aso_raw)
    pur_perc, longest_purine_stretch_raw = calculate_purine_stats(aso_raw)

    aso_name =  record.id

    added_row = pd.Series([aso_name, aso_raw.__str__(), maxcomp_raw, pur_perc, longest_purine_stretch_raw, None, None, None, None,
                            None, None, None, None],
                          index=output_df.columns)
    # do this but using contat and not append: output_df = output_df.append(added_row, ignore_index=True)
    output_df = pd.concat([output_df, added_row.to_frame().transpose()], ignore_index=True)

    seqs.append(SeqIO.SeqRecord(aso_target, aso_name, description=""))
    aso_seqs.append(SeqIO.SeqRecord(aso_raw, aso_name, description=""))


SeqIO.write(seqs, res_path + "/reference_sequences/aso_targets.fasta", "fasta")
SeqIO.write(aso_seqs, res_path + "/reference_sequences/shuffled_sequences.fasta", "fasta")
output_df.to_csv(res_path + "/outputs/result_table.tsv", sep="\t", index=False)

print("finished modift_PNAs script!")
