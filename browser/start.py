"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This will start the calculation in background
"""

import subprocess


def run_pipeline(script_path, args, done_file_path):
    """Generic pipeline runner: executes a shell script and writes a sentinel file on completion."""
    subprocess.run(["bash", script_path] + args)
    open(done_file_path, "w+").close()


def start_calculation(mason, fasta, gff, targets, length,  b_before, id, result_id, screen):
    """Start the MASON calculation."""
    run_pipeline(mason,
                 ["-f", fasta, "-g", gff, "-t", targets, "-l", length,
                  "-i", id, "-b", b_before, "-s", screen],
                 "./pnag/static/data/" + result_id + "/../done.txt")
    print("finished mason calculation!")


def start_scrambler(scrambler, fasta, gff, id, result_id, pna_input):
    """Start the scrambler calculation."""
    run_pipeline(scrambler,
                 ["-f", fasta, "-g", gff, "-i", id, "-p", pna_input],
                 "./pnag/static/data/" + result_id + "/done.txt")
    print("finished scrambler calculation!")


def start_checker(checker, fasta, gff, id, result_id, pna_input, checker_input):
    """Start the ASO-Checker calculation."""
    run_pipeline(checker,
                 ["-f", fasta, "-g", gff, "-i", id, "-p", pna_input, "-m", checker_input],
                 "./pnag/static/data/" + result_id + "/done.txt")
    print("finished ASO-Checker calculation!")
