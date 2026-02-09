import sys
from Bio import SeqIO
from cdifflib import CSequenceMatcher
import pandas as pd
from pna_utils import calculate_purine_stats, calculate_self_complementarity

seq_path = sys.argv[1]
res_path = sys.argv[2]
print(seq_path, "HIHIHIIHIH")

seqs = []
output_df = pd.DataFrame(columns=["ASO", "ASO_seq", "SC_bases", "pur_perc", "long_pur_stretch", "OT_tot", "OT_TIR"])

for record in SeqIO.parse(seq_path, "fasta"):
    print(record.id)
    print(record.seq)
    aso = record.seq

    maxcomp, aso_target = calculate_self_complementarity(aso)
    pur_perc, longest_purine_stretch = calculate_purine_stats(aso)

    aso_name = record.id

    added_row = pd.Series([aso_name, aso.__str__(), maxcomp, pur_perc, longest_purine_stretch, None, None],
                          index=output_df.columns)
    output_df = output_df.append(added_row, ignore_index=True)

    print(maxcomp)
    seqs.append(SeqIO.SeqRecord(aso_target, record.id, description=""))


SeqIO.write(seqs, res_path + "/reference_sequences/aso_targets.fasta", "fasta")
output_df.to_csv(res_path + "/outputs/result_table.tsv", sep="\t", index=False)
