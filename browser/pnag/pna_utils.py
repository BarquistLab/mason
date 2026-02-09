"""Shared utility functions for PNA/ASO sequence analysis across all tools."""

from cdifflib import CSequenceMatcher


def calculate_purine_stats(sequence):
    """Calculate purine percentage and longest purine stretch for a sequence.

    Args:
        sequence: DNA sequence (str or Bio.Seq object)

    Returns:
        tuple: (purine_percentage_str, longest_purine_stretch)
            purine_percentage_str is formatted as "XX.XX"
    """
    seq_str = str(sequence)
    pur = 0
    longest_purine_stretch = 0
    curstretch = 0
    for base in seq_str:
        if base in ("A", "G"):
            curstretch += 1
            pur += 1
            if curstretch > longest_purine_stretch:
                longest_purine_stretch += 1
        else:
            curstretch = 0
    pur_perc = "{:.2f}".format((pur / len(seq_str)) * 100)
    return pur_perc, longest_purine_stretch


def calculate_self_complementarity(aso_seq):
    """Calculate self-complementarity (longest match between ASO and its reverse complement).

    Args:
        aso_seq: Bio.Seq object of the ASO sequence

    Returns:
        tuple: (maxcomp, aso_target) where maxcomp is the longest match size
            and aso_target is the reverse complement Bio.Seq
    """
    aso_target = aso_seq.reverse_complement()
    maxcomp = CSequenceMatcher(None, aso_seq, aso_target).find_longest_match(
        0, len(aso_seq), 0, len(aso_target)).size
    return maxcomp, aso_target
