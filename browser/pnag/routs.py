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
from pnag.forms import startForm, startautoForm, ScrambledForm, CheckerForm
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


@app.route("/")
@app.route("/home")
def home():
    return render_template("home.html", title="Home")


@app.route("/start", methods=['GET', 'POST'])
def start():
    # form to upload the input files
    form = startForm()
    choices = [('upload', 'Own files')]
    choices += [(key, PRESETS[key][2]) for key in PRESETS.keys()]

    form.presets.choices = choices
    if form.validate_on_submit():
        paths = {}
        time_string = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S')

        additional_screen = request.form['add_screen']
        target_genes = form.genes.data
        target_genes = [x.strip() for x in target_genes.split(',')]
        print(target_genes)
        if "" in target_genes:
            target_genes.remove("")
        b_before = form.bases_before.data

        if len(target_genes) > 1:
            if form.essential.data is not None:
                target_genes += [form.essential.data]  # add space only if a lt exists:
        else:
            if form.essential.data is not None:
                target_genes += [form.essential.data]
        # create directory for result:
        os.mkdir("./pnag/static/data/"+time_string)
        # save uploaded data with timestring attached:
        if form.presets.data == 'upload':
            genome_file = os.path.join(app.root_path, 'static/data/', time_string + "/", form.genome.data.filename)
            gff_file = os.path.join(app.root_path, 'static/data/', time_string + "/", form.gff.data.filename)
            form.genome.data.save(genome_file)
            form.gff.data.save(gff_file)
            paths['genome'] = genome_file
            paths['gff'] = gff_file
        else:
            # genome first value in list
            files = PRESETS[form.presets.data]
            paths['genome'] = files[0]
            paths['gff'] = files[1]
        with open("./pnag/static/data/"+time_string+"/inputs.txt", "w+") as input_file:
            input_file.write(paths['genome'] + "," + paths['gff'])
        result_custom_id = form.custom_id.data

        # add selfcomp_errorfile
        efilename = './pnag/static/data/' +time_string+ "/error.txt"
        open(efilename, "w").close()

        # Now run MASON as background process while continuing with start.html and showing the "waiting" html:
        print(target_genes)
        for tgene in target_genes:
            print(tgene)
            resultid = time_string + "/" + tgene
            threading.Thread(target=start_calculation, name="masons", args=[path_parent.__str__() + "/mason.sh", paths['genome'],
                                                                            paths['gff'], tgene, str(form.len_PNA.data),
                                                                            str(form.mismatches.data), str(b_before),
                                                                            resultid, resultid, additional_screen]).start()
        with open("./pnag/static/data/"+time_string + "/inputs.txt", "a") as inputfile:
            inputfile.write("\n" + result_custom_id)
            inputfile.write("\n" + "; ".join(target_genes))
            inputfile.write("\n" + str(form.mismatches.data))
            inputfile.write("\n" + additional_screen)

        return redirect(url_for('result', result_id=time_string))
    return render_template("start.html", title="Start", essential=ESSENTIAL_GENES, form=form)


@app.route("/startauto", methods=['GET', 'POST'])
def startauto():
    # form to upload the input files
    form = startautoForm()
    choices = []
    choices += [(key, PRESETS[key][2]) for key in PRESETS.keys()]
    choices += [('upload', 'Own files')]

    form.presets.choices = choices
    if form.validate_on_submit():
        paths = {}
        time_string = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S')

        print(time_string)
        additional_screen = request.form['add_screen']
        target_genes = form.genes.data
        target_genes = [x.strip() for x in target_genes.split(',')]
        if "" in target_genes:
            target_genes.remove("")
        print(target_genes)
        b_before = form.bases_before.data

        if len(target_genes) > 1:
            if form.essential.data is not None:
                target_genes += [form.essential.data]  # add space only if a lt exists:
        else:
            if form.essential.data is not None:
                target_genes += [form.essential.data]
        # create directory for result:
        os.mkdir("./pnag/static/data/"+time_string)
        # save uploaded data with timestring attached:
        if form.presets.data == 'upload':
            genome_file = os.path.join(app.root_path, 'static/data/', time_string + "/", form.genome.data.filename)
            gff_file = os.path.join(app.root_path, 'static/data/', time_string + "/", form.gff.data.filename)
            form.genome.data.save(genome_file)
            form.gff.data.save(gff_file)
            paths['genome'] = genome_file
            paths['gff'] = gff_file
        else:
            # genome first value in list
            files = PRESETS[form.presets.data]
            paths['genome'] = files[0]
            paths['gff'] = files[1]
        with open("./pnag/static/data/"+time_string+"/inputs.txt", "w+") as input_file:
            input_file.write(paths['genome'] + "," + paths['gff'])
        result_custom_id = form.custom_id.data

        # add selfcomp_errorfile
        efilename = './pnag/static/data/' +time_string+ "/error.txt"
        open(efilename, "w").close()

        # Now run MASON as background process while continuing with start.html and showing the "waiting" html:
        for tgene in target_genes:
            resultid = time_string + "/" + tgene
            threading.Thread(target=start_calculation, name="masons", args=[path_parent.__str__() + "/mason.sh", paths['genome'],
                                                                            paths['gff'], tgene, str(form.len_PNA.data),
                                                                            str(form.mismatches.data), str(b_before),
                                                                            resultid, resultid, additional_screen]).start()
        with open("./pnag/static/data/"+time_string + "/inputs.txt", "a") as inputfile:
            inputfile.write("\n" + result_custom_id)
            inputfile.write("\n" + "; ".join(target_genes))
            inputfile.write("\n" + str(form.mismatches.data))
            inputfile.write("\n" + additional_screen)

        return redirect(url_for('result', result_id=time_string))
    return render_template("startauto.html", title="Start", essential=ESSENTIAL_GENES, form=form)


@app.route("/result/<result_id>")
def result(result_id):
    # each result gets its own page to access it. Just by calling .../result/<id of result>.
    time = re.sub("([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)",
                  "Date: \\1-\\2-\\3 ; Time: \\4:\\5:\\6", result_id)
    # get download paths:
    f = open('./pnag/static/data/' + result_id + "/inputs.txt")
    lines = f.readlines()
    custom_id = lines[1]
    tgenes = lines[2]
    mismatches = lines[3]
    add_screen = lines[4]
    genome_file, gff_file = lines[0].split(sep=",")
    f.close()

    dir_out = "../static/data/" + result_id
    rfinished = os.path.isfile("./pnag/static/data/" + result_id + "/done.txt")
    errs = "no errors"

    efilename = './pnag/static/data/' + result_id + "/error.txt"
    ltags_noPNAs = open(efilename).read().splitlines()
    errs = ltags_noPNAs

    all_output_dirs = []
    for dirs in os.listdir("./pnag/static/data/" + result_id):
        if "." not in dirs:
            all_output_dirs += [dirs]
    print(all_output_dirs)
    print("../" + "/".join(genome_file.split(sep="/")[-4:]))
    ffile = "../" + "/".join(genome_file.split(sep="/")[-4:])
    gfffile = "../" + "/".join(gff_file.split(sep="/")[-4:])
    print(all_output_dirs)
    # show results page
    return render_template("result.html", title="Result", result=result_id, dir_out=dir_out,
                           all_output_dirs=all_output_dirs,
                           rfin=rfinished, genome_file=ffile, gff_file=gfffile, custom_id=custom_id,
                           time=time, tgenes=tgenes, mismatches = mismatches, add_screen=add_screen,
                           errs=errs)


@app.route("/result_scrambler/<result_id>")
def result_scrambler(result_id):
    # each result gets its own page to access it. Just by calling .../result/<id of result>.
    time = re.sub("([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)",
                  "Date: \\1-\\2-\\3 | Time: \\4:\\5:\\6 h", result_id)
    # get download paths:
    f = open('./pnag/static/data/' + result_id + "/inputs.txt")
    lines = f.readlines()
    custom_id = lines[1]
    pnaseq = lines[2]
    genome_file, gff_file = lines[0].split(sep=",")
    f.close()

    dir_out = "../static/data/" + result_id
    rfinished = os.path.isfile("./pnag/static/data/" + result_id + "/done.txt")

    efilename = './pnag/static/data/' + result_id + "/error.txt"
    ltags_noPNAs = open(efilename).read().splitlines()
    errs = ltags_noPNAs

    all_output_dirs = []
    for dirs in os.listdir("./pnag/static/data/" + result_id):
        if "." not in dirs:
            all_output_dirs += [dirs]
    ffile = "../" + "/".join(genome_file.split(sep="/")[-4:])
    gfffile = "../" + "/".join(gff_file.split(sep="/")[-4:])
    # show results page
    return render_template("scrambler_result.html", title="Result", result=result_id, dir_out=dir_out,
                           all_output_dirs=all_output_dirs,
                           rfin=rfinished, genome_file=ffile, gff_file=gfffile, custom_id=custom_id, pnaseq=pnaseq,
                           time=time, errs=errs)


@app.route("/result_checker/<result_id>")
def result_checker(result_id):
    # each result gets its own page to access it. Just by calling .../result/<id of result>.
    time = re.sub("([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)",
                  "Date: \\1-\\2-\\3 | Time: \\4:\\5:\\6 h", result_id)
    # get download paths:
    f = open('./pnag/static/data/' + result_id + "/inputs.txt")
    lines = f.readlines()
    custom_id = lines[1]
    genome_file, gff_file = lines[0].split(sep=",")
    # get fasta pna sequence (make 1 string from list)
    pnaseq = lines[2:]
    pnaseq = "".join(pnaseq)
    f.close()

    dir_out = "../static/data/" + result_id
    rfinished = os.path.isfile("./pnag/static/data/" + result_id + "/done.txt")

    efilename = './pnag/static/data/' + result_id + "/error.txt"
    ltags_noPNAs = open(efilename).read().splitlines()
    errs = ltags_noPNAs

    all_output_dirs = []
    for dirs in os.listdir("./pnag/static/data/" + result_id):
        if "." not in dirs:
            all_output_dirs += [dirs]
    ffile = "../" + "/".join(genome_file.split(sep="/")[-4:])
    gfffile = "../" + "/".join(gff_file.split(sep="/")[-4:])
    # show results page
    return render_template("checker_result.html", title="Result ASO-Checker", result=result_id, dir_out=dir_out,
                           all_output_dirs=all_output_dirs,
                           rfin=rfinished, genome_file=ffile, gff_file=gfffile, custom_id=custom_id, pnaseq=pnaseq,
                           time=time, errs=errs)


@app.route("/delete_result/<result_id>")
def delete_result(result_id):
    shutil.rmtree("./pnag/static/data/" + result_id)
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


@app.route("/scrambler", methods=['GET', 'POST'])
def scrambler():
    form = ScrambledForm()
    print("abc")
    choices = [('upload', 'Own files')]
    choices += [(key, PRESETS[key][2]) for key in PRESETS.keys()]
    print(choices)
    form.presets.choices = choices
    form.seq_input
    if form.validate_on_submit():
        paths = {}
        time_string = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S')
        print(time_string)
        # create directory for result:
        os.mkdir("./pnag/static/data/" + time_string)
        # save uploaded data with timestring attached:
        if form.presets.data == 'upload':
            genome_file = os.path.join(app.root_path, 'static/data/', time_string + "/", form.genome.data.filename)
            gff_file = os.path.join(app.root_path, 'static/data/', time_string + "/", form.gff.data.filename)
            form.genome.data.save(genome_file)
            form.gff.data.save(gff_file)
            paths['genome'] = genome_file
            paths['gff'] = gff_file
        else:
            # genome first value in list
            files = PRESETS[form.presets.data]
            paths['genome'] = files[0]
            paths['gff'] = files[1]
        with open("./pnag/static/data/" + time_string + "/inputs.txt", "w+") as input_file:
            input_file.write(paths['genome'] + "," + paths['gff'])
        result_custom_id = form.custom_id.data
        pna_seq = form.seq_input.data
        pna_file_string = "./pnag/static/data/"+time_string+ "/pna_input.fasta"

        # add input PNA file:
        print("PNA sequence:")
        print(pna_seq)
        with open(pna_file_string, "w") as pna_file:
            pna_file.write(">PNA\n" + pna_seq)


        # add selfcomp_errorfile
        efilename = './pnag/static/data/' + time_string + "/error.txt"
        open(efilename, "w").close()

        # Now run scrambler as background process while continuing with scrambler.html and showing the "waiting" html:
        threading.Thread(target=start_scrambler, name="scramblers",
                         args=[path_parent.__str__() + "/scrambler.sh", paths['genome'],
                               paths['gff'], time_string, time_string, pna_seq]).start()
        with open("./pnag/static/data/" + time_string + "/inputs.txt", "a") as inputfile:
            inputfile.write("\n" + result_custom_id + "\n" + pna_seq)

        return redirect(url_for('result_scrambler', result_id=time_string))
    return render_template("scrambler.html", title="Scrambler", form=form)


@app.route("/ASO_checker", methods=['GET', 'POST'])
def checker():
    form = CheckerForm()
    choices = [('upload', 'Own files')]
    choices += [(key, PRESETS[key][2]) for key in PRESETS.keys()]
    print(choices)
    print("seqinput:")
    print(form.seq_input.data)
    print("done")
    form.presets.choices = choices
    if form.validate_on_submit():
        paths = {}
        time_string = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S')
        print(time_string)
        # create directory for result:
        os.mkdir("./pnag/static/data/" + time_string)
        # save uploaded data with timestring attached:
        if form.presets.data == 'upload':
            genome_file = os.path.join(app.root_path, 'static/data/', time_string + "/", form.genome.data.filename)
            gff_file = os.path.join(app.root_path, 'static/data/', time_string + "/", form.gff.data.filename)
            form.genome.data.save(genome_file)
            form.gff.data.save(gff_file)
            paths['genome'] = genome_file
            paths['gff'] = gff_file
        else:
            # genome first value in list
            files = PRESETS[form.presets.data]
            paths['genome'] = files[0]
            paths['gff'] = files[1]
        with open("./pnag/static/data/" + time_string + "/inputs.txt", "w+") as input_file:
            input_file.write(paths['genome'] + "," + paths['gff'])
        result_custom_id = form.custom_id.data
        pna_seq = form.seq_input.data
        pna_file_string = "./pnag/static/data/"+time_string+ "/pna_input.fasta"

        # add input PNA file:
        print("PNA sequence:")
        print(pna_seq)
        with open(pna_file_string, "w") as pna_file:
            pna_file.write(">PNA\n" + pna_seq)


        # add selfcomp_errorfile
        efilename = './pnag/static/data/' + time_string + "/error.txt"
        open(efilename, "w").close()

        # Now run checker as background process while continuing with checker.html and showing the "waiting" html:
        threading.Thread(target=start_checker, name="checkers",
                         args=[path_parent.__str__() + "/scrambler.sh", paths['genome'],
                               paths['gff'], time_string, time_string, pna_seq, "checker"]).start()
        with open("./pnag/static/data/" + time_string + "/inputs.txt", "a") as inputfile:
            inputfile.write("\n" + result_custom_id + "\n" + pna_seq)

        return redirect(url_for('result_checker', result_id=time_string))
    return render_template("checker.html", title="ASO-Checker", form=form)

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
