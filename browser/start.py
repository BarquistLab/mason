"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This will start the calculation in background
"""

import subprocess
from pnag import db
from pnag.models import Result


def start_calculation(mason, fasta, gff, targets, length, mismatches, id, result_id):
    """ This function starts the calculation of MASON
    """
    subprocess.run(["sh", mason, "-f", fasta, "-g", gff,
                    "-t", targets, "-l", length, "-m", mismatches,
                    "-i", id])
    r = Result.query.get(result_id)
    r.finish = True
    db.session.commit()
    print("finished mason calculation!")



