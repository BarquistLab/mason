"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This defines the structure of the website with its routs and functions within.
It tells which template to use and validates submitted forms.
"""
import os
from datetime import datetime
import re
import shutil
import pandas as pd
from flask import render_template, url_for, flash, redirect, request, send_file
# import needed things from other files in this package
from pnag import app, bcrypt, mail
from pnag.forms import startForm, ScrambledForm, CheckerForm
from flask_login import current_user
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
    return [('upload', 'Own files')] + [(key, PRESETS[key][2]) for key in PRESETS.keys()]


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

        if form.essential.data:
            target_genes += [form.essential.data]

        paths = _resolve_genome_paths(form, time_string)
        base_dir = _init_run_directory(time_string, paths)

        result_custom_id = form.custom_id.data

        for tgene in target_genes:
            resultid = time_string + "/" + tgene
            threading.Thread(target=start_calculation, name="masons",
                             args=[str(path_parent) + "/mason.sh", paths['genome'],
                                   paths['gff'], tgene, str(form.len_PNA.data),
                                   str(b_before),
                                   resultid, resultid, additional_screen]).start()

        with open(os.path.join(base_dir, "inputs.txt"), "a") as inputfile:
            inputfile.write("\n" + result_custom_id)
            inputfile.write("\n" + "; ".join(target_genes))
            inputfile.write("\n" + additional_screen)

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
        paths = _resolve_genome_paths(form, time_string)
        base_dir = _init_run_directory(time_string, paths)

        result_custom_id = form.custom_id.data
        pna_seq = form.seq_input.data
        pna_file_string = os.path.join(base_dir, "pna_input.fasta")

        with open(pna_file_string, "w") as pna_file:
            pna_file.write(">PNA\n" + pna_seq)

        threading.Thread(target=start_scrambler, name="scramblers",
                         args=[str(path_parent) + "/scrambler.sh", paths['genome'],
                               paths['gff'], time_string, time_string, pna_seq]).start()

        with open(os.path.join(base_dir, "inputs.txt"), "a") as inputfile:
            inputfile.write("\n" + result_custom_id + "\n" + pna_seq)

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
        paths = _resolve_genome_paths(form, time_string)
        base_dir = _init_run_directory(time_string, paths)

        result_custom_id = form.custom_id.data
        pna_seq = form.seq_input.data
        pna_file_string = os.path.join(base_dir, "pna_input.fasta")

        with open(pna_file_string, "w") as pna_file:
            pna_file.write(">PNA\n" + pna_seq)

        threading.Thread(target=start_checker, name="checkers",
                         args=[str(path_parent) + "/scrambler.sh", paths['genome'],
                               paths['gff'], time_string, time_string, pna_seq, "checker"]).start()

        with open(os.path.join(base_dir, "inputs.txt"), "a") as inputfile:
            inputfile.write("\n" + result_custom_id + "\n" + pna_seq)

        return redirect(url_for('result_checker', result_id=time_string))
    return render_template("checker.html", title="ASO-Checker", form=form)


@app.route("/result/<result_id>")
def result(result_id):
    ctx = _read_result_context(result_id)
    tgenes = ctx['lines'][2]
    add_screen = ctx['lines'][3]

    return render_template("result.html", title="Result", result=result_id, dir_out=ctx['dir_out'],
                           all_output_dirs=ctx['all_output_dirs'],
                           rfin=ctx['rfinished'], genome_file=ctx['ffile'], gff_file=ctx['gfffile'],
                           custom_id=ctx['custom_id'],
                           time=ctx['time'], tgenes=tgenes, add_screen=add_screen,
                           errs=ctx['errs'])


@app.route("/result_scrambler/<result_id>")
def result_scrambler(result_id):
    ctx = _read_result_context(result_id)
    pnaseq = ctx['lines'][2]

    return render_template("scrambler_result.html", title="Result", result=result_id, dir_out=ctx['dir_out'],
                           all_output_dirs=ctx['all_output_dirs'],
                           rfin=ctx['rfinished'], genome_file=ctx['ffile'], gff_file=ctx['gfffile'],
                           custom_id=ctx['custom_id'], pnaseq=pnaseq,
                           time=ctx['time'], errs=ctx['errs'])


@app.route("/result_checker/<result_id>")
def result_checker(result_id):
    ctx = _read_result_context(result_id)
    # get fasta pna sequence (make 1 string from list)
    pnaseq = "".join(ctx['lines'][2:])

    return render_template("checker_result.html", title="Result ASO-Checker", result=result_id,
                           dir_out=ctx['dir_out'],
                           all_output_dirs=ctx['all_output_dirs'],
                           rfin=ctx['rfinished'], genome_file=ctx['ffile'], gff_file=ctx['gfffile'],
                           custom_id=ctx['custom_id'], pnaseq=pnaseq,
                           time=ctx['time'], errs=ctx['errs'])


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
    first.split('_')
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


def save_file(form_file, prefix, time_str):
    # saves files with prefix and time in data folder
    if form_file:
        if prefix == 'genes':
            file_fn = prefix + time_str + str(current_user.id) + '.txt'
            file_path = os.path.join(app.root_path, 'static/data', file_fn)
            f = open(file_path, 'w')
            f.write(form_file)
            f.close()
        else:
            _, f_ext = os.path.splitext(form_file.filename)
            file_fn = prefix + time_str + str(current_user.id) + f_ext
            file_path = os.path.join(app.root_path, 'static/data', file_fn)
            form_file.save(file_path)
    elif prefix == 'pnas':
        file_fn = prefix + time_str + str(current_user.id) + '.fasta'
        file_path = os.path.join(app.root_path, 'static/data', file_fn)
    else:
        file_fn = prefix + time_str + str(current_user.id) + '.fasta'
        file_path = os.path.join(path_parent, 'pnag/static/data/', file_fn)
    return file_path


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
