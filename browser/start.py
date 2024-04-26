"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This will start the calculation in background
"""

import subprocess


def start_calculation(mason, fasta, gff, targets, length, mismatches, b_before,  id, result_id, screen):
    """ This function starts the calculation of MASON
    """
    print(b_before)
    print(result_id)
    print(id)
        
    subprocess.run(["bash", mason, "-f", fasta, "-g", gff,
                    "-t", targets, "-l", length, "-m", mismatches,
                    "-i", id, "-b", b_before, "-s", screen])
    f = open("./pnag/static/data/" + result_id + "/../done.txt", "w+")
    f.close()
    print("finished mason calculation!")


def start_scrambler(scrambler, fasta, gff, id, result_id, pna_input):
    """ This function starts the calculation of scrambler
    """

    subprocess.run(["bash", scrambler, "-f", fasta, "-g", gff, "-i", id, "-p", pna_input])
    print(result_id)
    f = open("./pnag/static/data/" + result_id + "/done.txt", "w+")
    f.close()
    print("finished scrambler calculation!")


def start_checker(checker, fasta, gff, id, result_id, pna_input, checker_input):
    """ This function starts the calculation of scrambler
    """

    subprocess.run(["bash", checker, "-f", fasta, "-g", gff, "-i", id, "-p", pna_input, "-m", checker_input])
    print(result_id)
    f = open("./pnag/static/data/" + result_id + "/done.txt", "w+")
    f.close()
    print("finished ASO-Checker calculation!")


