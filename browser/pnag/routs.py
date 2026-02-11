"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This defines the structure of the website with its routs and functions within.
It tells which template to use and validates submitted forms.
"""
import os
import glob
import subprocess
import zipfile
from datetime import datetime
import re
import shutil
from flask import render_template, url_for, flash, redirect, request, send_file
from pnag import app
from pnag.forms import startForm, ScrambledForm, CheckerForm
import base64
import json
import threading
from start import start_calculation, start_scrambler, start_checker
from pathlib import Path

path_parent = Path(app.root_path).parent

p = os.path.join(app.root_path, 'static/data/presets/')
PRESETS = {
    'e_coli': [p + 'e_coli_K12.fna', p + 'e_coli_K12.gff', 'E. coli str. K-12 substr. MG1655'],
    's_typhi': [p + 's_typhi.fa', p + 's_typhi.gff', 'Salmonella enterica subsp. enterica serovar Typhimurium SL1344'],
    'c_diffi': [p + 'c_diffi630.fasta', p + 'c_diffi630.gff3', 'Clostridium difficile 630'],
    'Fuso': [p + 'Fuso.fasta', p + 'Fuso.gff3', 'Fusobacterium nucleatum ATCC 23726']
}
with open('ESSENTIAL_GENES.json', 'r') as f:
    ESSENTIAL_GENES = json.load(f)


def _preset_choices():
    """Build the preset choices list used by all forms."""
    return ([('upload', 'Own files'), ('ncbi', 'NCBI assembly accession')]
            + [(key, PRESETS[key][2]) for key in PRESETS.keys()])


def _download_ncbi_accession(accession, dest_dir):
    """Download genome FASTA and GFF3 from NCBI using the datasets CLI.

    Returns dict with 'genome' and 'gff' paths.
    Raises RuntimeError on failure.
    """
    os.makedirs(dest_dir, exist_ok=True)
    zip_path = os.path.join(dest_dir, 'ncbi_dataset.zip')

    result = subprocess.run(
        ['datasets', 'download', 'genome', 'accession', accession,
         '--include', 'gff3,genome', '--filename', zip_path, '--no-progressbar'],
        capture_output=True, text=True, timeout=300
    )
    if result.returncode != 0:
        raise RuntimeError(
            f'Failed to download accession {accession}. '
            'Please check that the accession is valid (e.g. GCF_000005845.2).'
        )

    try:
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(dest_dir)
    except zipfile.BadZipFile:
        raise RuntimeError(f'Downloaded file for {accession} is not a valid zip archive.')
    finally:
        if os.path.exists(zip_path):
            os.remove(zip_path)

    data_dir = os.path.join(dest_dir, 'ncbi_dataset', 'data', accession)
    fna_files = glob.glob(os.path.join(data_dir, '*.fna'))
    gff_files = glob.glob(os.path.join(data_dir, '*.gff'))

    if not fna_files or not gff_files:
        raise RuntimeError(
            f'Could not find FASTA/GFF files for accession {accession}. '
            'The assembly may not include both genome and annotation files.'
        )

    return {'genome': fna_files[0], 'gff': gff_files[0]}


def _resolve_genome_paths(form, time_string):
    """Handle preset-or-upload file logic. Returns dict with 'genome' and 'gff' paths."""
    paths = {}
    if form.presets.data == 'upload':
        # Create directory first using proper path joining
        base_dir = os.path.join(app.root_path, 'static/data', time_string)
        os.makedirs(base_dir, exist_ok=True)

        # Now save files to the created directory
        genome_file = os.path.join(base_dir, form.genome.data.filename)
        gff_file = os.path.join(base_dir, form.gff.data.filename)

        form.genome.data.save(genome_file)
        form.gff.data.save(gff_file)
        paths['genome'] = genome_file
        paths['gff'] = gff_file
    elif form.presets.data == 'ncbi':
        base_dir = os.path.join(app.root_path, 'static/data', time_string)
        paths = _download_ncbi_accession(form.ncbi_accession.data.strip(), base_dir)
    else:
        files = PRESETS[form.presets.data]
        paths['genome'] = files[0]
        paths['gff'] = files[1]
    return paths


def _init_run_directory(time_string, paths):
    """Create timestamped result directory, write inputs.txt and error.txt."""
    base_dir = os.path.join(app.root_path, 'static/data', time_string)
    os.makedirs(base_dir, exist_ok=True)
    with open(os.path.join(base_dir, "inputs.txt"), "w+") as input_file:
        input_file.write(paths['genome'] + "," + paths['gff'])
    open(os.path.join(base_dir, "error.txt"), "w").close()
    return base_dir


def _format_checker_sequences(seq_input):
    """Convert plain sequences to FASTA format with aso_NNN headers if needed.
    Also uppercases all sequence characters."""
    lines = [l.strip() for l in seq_input.strip().splitlines() if l.strip()]
    if not lines:
        return seq_input
    if lines[0].startswith('>'):
        # Already FASTA — uppercase sequence lines only
        return '\n'.join(l if l.startswith('>') else l.upper()
                         for l in seq_input.strip().splitlines())
    # Plain sequences — wrap each line with an aso_NNN header
    fasta_lines = []
    for i, line in enumerate(lines, 1):
        fasta_lines.append(f'>aso_{i:03d}')
        fasta_lines.append(line.upper())
    return '\n'.join(fasta_lines)


def _read_result_context(result_id):
    """Read shared result page context: inputs.txt, done status, output dirs, file paths."""
    base_dir = os.path.join(app.root_path, 'static/data', result_id)
    time = re.sub("([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)",
                  "Date: \\1-\\2-\\3 | Time: \\4:\\5:\\6 h", result_id)

    with open(os.path.join(base_dir, "inputs.txt")) as f:
        lines = f.readlines()

    genome_file, gff_file = lines[0].split(sep=",")
    custom_id = lines[1]

    dir_out = "../static/data/" + result_id
    rfinished = os.path.isfile(os.path.join(base_dir, "done.txt"))

    errs = open(os.path.join(base_dir, "error.txt")).read().splitlines()

    all_output_dirs = []
    for dirs in os.listdir(base_dir):
        if "." not in dirs:
            all_output_dirs += [dirs]

    ffile = "../" + "/".join(genome_file.split(sep="/")[-4:])
    gfffile = "../" + "/".join(gff_file.split(sep="/")[-4:])

    return {
        'lines': lines,
        'time': time,
        'custom_id': custom_id,
        'dir_out': dir_out,
        'rfinished': rfinished,
        'errs': errs,
        'all_output_dirs': all_output_dirs,
        'ffile': ffile,
        'gfffile': gfffile,
    }


@app.route("/")
@app.route("/home")
def home():
    return render_template("home.html", title="Home")


@app.route("/startauto", methods=['GET', 'POST'])
@app.route("/start", methods=['GET', 'POST'])
def start():
    form = startForm()
    form.presets.choices = _preset_choices()

    # Autofill only on initial GET for the auto route
    if request.method == 'GET' and request.path == '/startauto':
        form.custom_id.data = "example_auto_001"
        form.presets.data = "e_coli"
        form.genes.data = "b0081"
        form.len_PNA.data = 10
        form.bases_before.data = ""
        form.essential.data = []

    if form.validate_on_submit():
        time_string = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S')
        additional_screen = request.form['add_screen']
        target_genes = [x.strip() for x in form.genes.data.split(',') if x.strip()]
        b_before = form.bases_before.data
        use_ml = "yes" if form.use_ml.data else "no"

        if form.essential.data:
            target_genes += [form.essential.data]

        try:
            paths = _resolve_genome_paths(form, time_string)
        except RuntimeError as e:
            flash(str(e), 'danger')
            return render_template("start.html", title="Start", essential=ESSENTIAL_GENES, form=form)
        base_dir = _init_run_directory(time_string, paths)

        result_custom_id = form.custom_id.data

        open(os.path.join(path_parent, "logfile_masonscript.log"), "w").close()

        for tgene in target_genes:
            resultid = time_string + "/" + tgene
            threading.Thread(target=start_calculation, name="masons",
                             args=[str(path_parent) + "/mason.sh", paths['genome'],
                                   paths['gff'], tgene, str(form.len_PNA.data),
                                   str(b_before),
                                   resultid, resultid, additional_screen,
                                   use_ml]).start()

        with open(os.path.join(base_dir, "inputs.txt"), "a") as inputfile:
            inputfile.write("\n" + result_custom_id)
            inputfile.write("\n" + "; ".join(target_genes))
            inputfile.write("\n" + additional_screen)
            inputfile.write("\n" + use_ml)

        return redirect(url_for('result', result_id=time_string))
    return render_template("start.html", title="Start", essential=ESSENTIAL_GENES, form=form)


@app.route("/scramblerauto", methods=['GET', 'POST'])
@app.route("/scrambler", methods=['GET', 'POST'])
def scrambler():
    form = ScrambledForm()
    form.presets.choices = _preset_choices()

    # Autofill only on GET for the auto route
    if request.method == 'GET' and request.path == '/scramblerauto':
        form.custom_id.data = "example_scramble"
        form.presets.data = "e_coli"
        form.seq_input.data = "ATCTCGCAT"

    if form.validate_on_submit():
        time_string = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S')
        additional_screen = request.form['add_screen']
        try:
            paths = _resolve_genome_paths(form, time_string)
        except RuntimeError as e:
            flash(str(e), 'danger')
            return render_template("scrambler.html", title="Scrambler", form=form)
        base_dir = _init_run_directory(time_string, paths)

        result_custom_id = form.custom_id.data
        pna_seq = form.seq_input.data.upper()
        pna_file_string = os.path.join(base_dir, "pna_input.fasta")

        with open(pna_file_string, "w") as pna_file:
            pna_file.write(">PNA\n" + pna_seq)

        open(os.path.join(path_parent, "logfile_masonscript.log"), "w").close()

        threading.Thread(target=start_scrambler, name="scramblers",
                         args=[str(path_parent) + "/scrambler.sh", paths['genome'],
                               paths['gff'], time_string, time_string, pna_seq,
                               additional_screen]).start()

        with open(os.path.join(base_dir, "inputs.txt"), "a") as inputfile:
            inputfile.write("\n" + result_custom_id + "\n" + pna_seq)
            inputfile.write("\n" + additional_screen)

        return redirect(url_for('result_scrambler', result_id=time_string))
    return render_template("scrambler.html", title="Scrambler", form=form)


@app.route("/checkerauto", methods=['GET', 'POST'])
@app.route("/ASO_checker", methods=['GET', 'POST'])
def checker():
    form = CheckerForm()
    form.presets.choices = _preset_choices()

    # Autofill only on GET for the auto route
    if request.method == 'GET' and request.path == '/checkerauto':
        form.custom_id.data = "example_checker_001"
        form.presets.data = "e_coli"
        form.seq_input.data = """>ASO_1
ATGCTGCTGC
>ASO_2
CGATACGTGA
>ASO_3
ATATATATA"""

    if form.validate_on_submit():
        time_string = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S')
        additional_screen = request.form['add_screen']
        target_gene = form.target_gene.data.strip() if form.target_gene.data else ""
        use_ml = "yes" if (target_gene and form.use_ml.data) else "no"
        try:
            paths = _resolve_genome_paths(form, time_string)
        except RuntimeError as e:
            flash(str(e), 'danger')
            return render_template("checker.html", title="ASO-Checker", form=form)
        base_dir = _init_run_directory(time_string, paths)

        result_custom_id = form.custom_id.data
        pna_seq = _format_checker_sequences(form.seq_input.data)
        pna_file_string = os.path.join(base_dir, "pna_input.fasta")

        with open(pna_file_string, "w") as pna_file:
            pna_file.write(pna_seq)

        open(os.path.join(path_parent, "logfile_masonscript.log"), "w").close()

        threading.Thread(target=start_checker, name="checkers",
                         args=[str(path_parent) + "/scrambler.sh", paths['genome'],
                               paths['gff'], time_string, time_string, pna_seq, "checker",
                               additional_screen],
                         kwargs={"target_gene": target_gene, "use_ml": use_ml}).start()

        with open(os.path.join(base_dir, "inputs.txt"), "a") as inputfile:
            inputfile.write("\n" + result_custom_id + "\n" + pna_seq)
            inputfile.write("\n" + additional_screen)
            if target_gene:
                inputfile.write("\n" + target_gene)
                inputfile.write("\n" + use_ml)

        return redirect(url_for('result_checker', result_id=time_string))
    return render_template("checker.html", title="ASO-Checker", form=form)


@app.route("/result/<result_id>")
def result(result_id):
    ctx = _read_result_context(result_id)
    tgenes = ctx['lines'][2]
    add_screen = ctx['lines'][3]
    use_ml = ctx['lines'][4].strip() if len(ctx['lines']) > 4 else "no"

    return render_template("result.html", title="Result", result=result_id, dir_out=ctx['dir_out'],
                           all_output_dirs=ctx['all_output_dirs'],
                           rfin=ctx['rfinished'], genome_file=ctx['ffile'], gff_file=ctx['gfffile'],
                           custom_id=ctx['custom_id'],
                           time=ctx['time'], tgenes=tgenes, add_screen=add_screen,
                           errs=ctx['errs'], use_ml=use_ml)


@app.route("/result_scrambler/<result_id>")
def result_scrambler(result_id):
    ctx = _read_result_context(result_id)
    pnaseq = ctx['lines'][2]
    add_screen = ctx['lines'][3].strip() if len(ctx['lines']) > 3 else "none"

    return render_template("scrambler_result.html", title="Result", result=result_id, dir_out=ctx['dir_out'],
                           all_output_dirs=ctx['all_output_dirs'],
                           rfin=ctx['rfinished'], genome_file=ctx['ffile'], gff_file=ctx['gfffile'],
                           custom_id=ctx['custom_id'], pnaseq=pnaseq,
                           time=ctx['time'], errs=ctx['errs'], add_screen=add_screen)


@app.route("/result_checker/<result_id>")
def result_checker(result_id):
    ctx = _read_result_context(result_id)
    lines = ctx['lines']

    # New format: lines are genome,gff / custom_id / pnaseq... / screen / target_gene / use_ml
    # Detect new format: last line is "yes"/"no" (use_ml) and second-to-last is target_gene
    target_gene = ""
    use_ml = "no"
    last_line = lines[-1].strip() if lines else ""
    second_last = lines[-2].strip() if len(lines) >= 2 else ""

    if last_line in ("yes", "no") and second_last not in ("none", "human", "microbiome"):
        # New format with target_gene and use_ml appended
        use_ml = last_line
        target_gene = second_last
        add_screen = lines[-3].strip() if len(lines) >= 3 else "none"
        pnaseq = "".join(lines[2:-3])
    elif last_line in ("none", "human", "microbiome"):
        add_screen = last_line
        pnaseq = "".join(lines[2:-1])
    else:
        add_screen = "none"
        pnaseq = "".join(lines[2:])

    # Read varna_positions.tsv for VARNA display when target gene is set
    varna_asos = []
    warnings = []
    if target_gene:
        base_dir = os.path.join(app.root_path, 'static/data', result_id)
        varna_path = os.path.join(base_dir, "outputs", "varna_positions.tsv")
        if os.path.isfile(varna_path):
            with open(varna_path) as vf:
                for line in vf:
                    parts = line.strip().split("\t")
                    if parts and parts[0] != "aso_name":
                        varna_asos.append(parts[0])
        warnings_path = os.path.join(base_dir, "outputs", "warnings.txt")
        if os.path.isfile(warnings_path):
            with open(warnings_path) as wf:
                warnings = [l.strip() for l in wf if l.strip()]

    return render_template("checker_result.html", title="Result ASO-Checker", result=result_id,
                           dir_out=ctx['dir_out'],
                           all_output_dirs=ctx['all_output_dirs'],
                           rfin=ctx['rfinished'], genome_file=ctx['ffile'], gff_file=ctx['gfffile'],
                           custom_id=ctx['custom_id'], pnaseq=pnaseq,
                           time=ctx['time'], errs=ctx['errs'], add_screen=add_screen,
                           target_gene=target_gene, use_ml=use_ml,
                           varna_asos=varna_asos, warnings=warnings)


@app.route("/delete_result/<result_id>")
def delete_result(result_id):
    shutil.rmtree(os.path.join(app.root_path, 'static/data', result_id))
    return redirect(url_for('home'))


@app.route("/download/<path>")
def download(path):
    # just for download the file in path
    path = url_decode(path)
    first, f_ext = os.path.splitext(path)
    preset = first
    preset = preset.split('/')
    if 'presets' in preset:
        return send_file(path, as_attachment=True)
    else:
        # just the owner can download the file
        try:
            return send_file(path, as_attachment=True)
        except FileNotFoundError:
            flash(u'This file does not exist. Maybe no target genes found to build PNA.fasta.', 'danger')
            return redirect(url_for('home'))


@app.route("/about")
def about():
    return render_template("about.html", title="About")


@app.route("/help")
def help():
    return render_template("help.html", title="Help")


@app.errorhandler(404)
def error_404(error):
    return render_template('errors/404.html'), 404


@app.errorhandler(403)
def error_403(error):
    return render_template('errors/403.html'), 403


@app.errorhandler(500)
def error_500(error):
    return render_template('errors/500.html'), 500


@app.template_filter()
def url_encode(env, string):
    urlSafeEncodedBytes = base64.urlsafe_b64encode(string.encode("utf-8"))
    urlSafeEncodedStr = str(urlSafeEncodedBytes, "utf-8")
    return urlSafeEncodedStr


def url_decode(code):
    decodedBytes = base64.urlsafe_b64decode(code)
    decodedStr = str(decodedBytes, "utf-8")
    return decodedStr
