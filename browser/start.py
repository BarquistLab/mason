"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This will start the calculation in background
"""

import subprocess
from pnag import db
from pnag.models import Result


def start_calculation(mason, fasta, gff, targets, length, mismatches, b_before,  id, result_id):
    """ This function starts the calculation of MASON
    """
    print(b_before)
        
    subprocess.run(["bash", mason, "-f", fasta, "-g", gff,
                    "-t", targets, "-l", length, "-m", mismatches,
                    "-i", id, "-b", b_before])
    r = Result.query.get(result_id)
    r.finish = True
    db.session.commit()
    print("finished mason calculation!")



