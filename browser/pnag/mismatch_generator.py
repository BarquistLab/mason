"""Generate mismatch control sequences by swapping pairs of non-equal bases.

Usage: python mismatch_generator.py <input_fasta> <res_path> <num_mismatches>

For N mismatches (must be even), generates all possible sequences obtained by
swapping N/2 non-overlapping pairs of positions where the bases differ.

Outputs the same files as scrambler_modify_pnas.py:
  - <res_path>/reference_sequences/aso_targets.fasta
  - <res_path>/reference_sequences/shuffled_sequences.fasta
  - <res_path>/outputs/result_table.tsv
"""
import sys
from itertools import combinations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from pna_utils import calculate_purine_stats, calculate_self_complementarity


def generate_mismatch_sequences(seq_str, num_mismatches):
    """Generate all mismatch sequences by swapping pairs of non-equal bases.

    Args:
        seq_str: The original sequence string.
        num_mismatches: Number of mismatches to introduce (must be even).

    Returns:
        List of (name, sequence_string) tuples.
    """
    seq_len = len(seq_str)
    num_swaps = num_mismatches // 2

    # Find all pairs of positions where bases differ
    diff_pairs = []
    for i in range(seq_len):
        for j in range(i + 1, seq_len):
            if seq_str[i] != seq_str[j]:
                diff_pairs.append((i, j))

    # Generate all combinations of num_swaps non-overlapping pairs
    results = []
    seen = set()

    for swap_combo in combinations(diff_pairs, num_swaps):
        # Check that all positions in this combination are unique (non-overlapping)
        positions = []
        for pair in swap_combo:
            positions.extend(pair)
        if len(set(positions)) != len(positions):
            continue

        # Apply the swaps
        seq_list = list(seq_str)
        for i, j in swap_combo:
            seq_list[i], seq_list[j] = seq_list[j], seq_list[i]

        new_seq = ''.join(seq_list)

        # Skip duplicates and the original sequence
        if new_seq in seen or new_seq == seq_str:
            continue

        # Filter: max consecutive matching bases with original must be < 6
        max_consec = 0
        consec = 0
        for a, b in zip(seq_str, new_seq):
            if a == b:
                consec += 1
                max_consec = max(max_consec, consec)
            else:
                consec = 0
        if max_consec >= 6:
            continue

        seen.add(new_seq)

        name = f"PNA_mm_{len(results) + 1:03d}"
        results.append((name, new_seq))

    return results


def main():
    input_fasta = sys.argv[1]
    res_path = sys.argv[2]
    num_mismatches = int(sys.argv[3])

    if num_mismatches % 2 != 0:
        print(f"Error: num_mismatches must be even, got {num_mismatches}", file=sys.stderr)
        sys.exit(1)

    # Read input sequence
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        print("Error: no sequences found in input FASTA", file=sys.stderr)
        sys.exit(1)

    original_seq = str(records[0].seq)
    print(f"Input sequence: {original_seq} (length {len(original_seq)})")
    print(f"Generating sequences with {num_mismatches} mismatches ({num_mismatches // 2} swaps)...")

    mismatch_seqs = generate_mismatch_sequences(original_seq, num_mismatches)
    print(f"Generated {len(mismatch_seqs)} unique mismatch sequences")

    # Build output dataframe and FASTA records (same structure as scrambler_modify_pnas.py)
    output_df = pd.DataFrame(columns=["ASO", "ASO_seq", "SC_bases", "long_pur_stretch",
                                       "OT_TIR_0mm", "OT_TIR_1mm", "OT_TIR_2mm", "OT_TIR_3mm",
                                       "OT_tot_0mm", "OT_tot_1mm", "OT_tot_2mm", "OT_tot_3mm"])

    aso_targets = []
    aso_seqs = []

    # First: add the input PNA
    aso_raw = records[0].seq
    maxcomp_raw, aso_target_raw = calculate_self_complementarity(aso_raw)
    _, longest_pur_raw = calculate_purine_stats(aso_raw)

    row = pd.Series(["input_PNA", str(aso_raw), maxcomp_raw, longest_pur_raw,
                      None, None, None, None, None, None, None, None],
                     index=output_df.columns)
    output_df = pd.concat([output_df, row.to_frame().transpose()], ignore_index=True)
    aso_targets.append(SeqRecord(aso_target_raw, id="input_PNA", description=""))
    aso_seqs.append(SeqRecord(aso_raw, id="input_PNA", description=""))

    # Then: add mismatch sequences
    for name, seq_str in mismatch_seqs:
        aso = Seq(seq_str)
        maxcomp, aso_target = calculate_self_complementarity(aso)
        _, longest_pur = calculate_purine_stats(aso)

        row = pd.Series([name, seq_str, maxcomp, longest_pur,
                          None, None, None, None, None, None, None, None],
                         index=output_df.columns)
        output_df = pd.concat([output_df, row.to_frame().transpose()], ignore_index=True)
        aso_targets.append(SeqRecord(aso_target, id=name, description=""))
        aso_seqs.append(SeqRecord(aso, id=name, description=""))

    # Write output files
    SeqIO.write(aso_targets, res_path + "/reference_sequences/aso_targets.fasta", "fasta")
    SeqIO.write(aso_seqs, res_path + "/reference_sequences/shuffled_sequences.fasta", "fasta")
    output_df.to_csv(res_path + "/outputs/result_table.tsv", sep="\t", index=False)

    print(f"Written {len(mismatch_seqs)} mismatch sequences + input PNA to output files")


if __name__ == "__main__":
    main()
