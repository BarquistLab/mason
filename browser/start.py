"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This will start the calculation in background
"""

import subprocess


def start_calculation(mason, fasta, gff, targets, length, mismatches, b_before,  id, result_id):
    """ This function starts the calculation of MASON
    """
    print(b_before)
    print(result_id)
    print(id)
        
    subprocess.run(["bash", mason, "-f", fasta, "-g", gff,
                    "-t", targets, "-l", length, "-m", mismatches,
                    "-i", id, "-b", b_before])
    f = open("./pnag/static/data/" + result_id + "/../done.txt", "w+")
    f.close()
    print("finished mason calculation!")



