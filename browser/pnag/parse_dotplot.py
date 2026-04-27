"""Parse ViennaRNA dot plot and emit per-base paired probabilities.

Usage: python parse_dotplot.py <dot.ps> <seq_length>

Reads the upper-triangle "ubox" lines (i j sqrt(p)) from a ViennaRNA dot plot
and writes, to stdout, semicolon-separated per-base paired probabilities
(P_paired[i] = sum_j P(i,j)) — formatted for VARNA's -colorMap argument.
"""
import re
import sys

dp_path = sys.argv[1]
seq_len = int(sys.argv[2])

paired = [0.0] * seq_len
ubox_re = re.compile(r'^\s*(\d+)\s+(\d+)\s+([\d.eE+-]+)\s+ubox\s*$')

with open(dp_path) as f:
    for line in f:
        m = ubox_re.match(line)
        if not m:
            continue
        i, j, sqrt_p = int(m.group(1)), int(m.group(2)), float(m.group(3))
        p = sqrt_p ** 2
        if 1 <= i <= seq_len:
            paired[i - 1] += p
        if 1 <= j <= seq_len:
            paired[j - 1] += p

paired = [min(1.0, max(0.0, p)) for p in paired]
print(";".join(f"{p:.4f}" for p in paired))
