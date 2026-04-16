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
import logging
from flask import render_template, url_for, flash, redirect, request, send_file
from pnag import app
from pnag.forms import startForm, ScrambledForm, CheckerForm
import base64
import json
import threading
from start import start_calculation, start_scrambler, start_checker
from pathlib import Path

# Persistent usage logger — one line per job submission
_usage_logger = logging.getLogger('mason_usage')
_usage_logger.setLevel(logging.INFO)
_usage_handler = logging.FileHandler(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static/data/usage.log'))
_usage_handler.setFormatter(logging.Formatter('%(message)s'))
_usage_logger.addHandler(_usage_handler)

path_parent = Path(app.root_path).parent

p = os.path.join(app.root_path, 'static/data/presets/')
PRESETS = {
    # Local preset genomes (bundled files)
    'e_coli':  {'genome': p + 'e_coli_K12.fna', 'gff': p + 'e_coli_K12.gff',
                'display': 'Escherichia coli str. K-12 substr. MG1655'},
    's_typhi': {'genome': p + 's_typhi.fa', 'gff': p + 's_typhi.gff',
                'display': 'Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344'},
    'c_diffi': {'genome': p + 'c_diffi630.fasta', 'gff': p + 'c_diffi630.gff3',
                'display': 'Clostridioides difficile str. 630'},
    'Fuso':    {'genome': p + 'Fuso.fasta', 'gff': p + 'Fuso.gff3',
                'display': 'Fusobacterium nucleatum subsp. nucleatum ATCC 23726'},
    # NCBI accession presets (genome downloaded on first use)
    # locus_tag_prefix: post-process GFF to use old_locus_tag matching this prefix
    # (NCBI RefSeq GFFs use re-annotated *_RS locus tags; originals are in old_locus_tag)
    'ecoli_bw25113':   {'accession': 'GCF_000750555.1',
                        'display': 'Escherichia coli str. BW25113',
                        'locus_tag_prefix': 'BW25113_'},
    'ecoli_ec958':     {'accession': 'GCF_000285655.3',
                        'display': 'Escherichia coli str. EC958 (ST131)',
                        'locus_tag_prefix': 'EC958_'},
    'ecoli_nctc13441': {'accession': 'GCF_900448475.1',
                        'display': 'Escherichia coli str. NCTC13441 (UPEC ST131)',
                        'locus_tag_prefix': 'NCTC13441_'},
    'c_rod_icc168':    {'accession': 'GCF_000027085.1',
                        'display': 'Citrobacter rodentium str. ICC168',
                        'locus_tag_prefix': 'ROD_'},
    'kpneu_ecl8':      {'accession': 'GCF_000315385.1',
                        'display': 'Klebsiella pneumoniae str. Ecl8',
                        'locus_tag_prefix': 'BN373_'},
    'kpneu_rh201207':  {'accession': 'GCA_905477585.1',
                        'display': 'Klebsiella pneumoniae str. RH201207',
                        'locus_tag_prefix': 'KPNRH_'},
    's_enteritidis':   {'accession': 'GCF_000009505.1',
                        'display': 'Salmonella enterica subsp. enterica serovar Enteritidis str. P125109',
                        'locus_tag_prefix': 'SEN'},
    's_typhi_ty2':     {'accession': 'GCF_000007545.1',
                        'display': 'Salmonella enterica subsp. enterica serovar Typhi str. Ty2',
                        'locus_tag_prefix': 't'},
    's_typh_a130':     {'accession': 'GCF_000027025.1',
                        'display': 'Salmonella enterica subsp. enterica serovar Typhimurium str. A130',
                        'locus_tag_prefix': 'STM_MW'},
    's_typh_d23580':   {'accession': 'GCF_000027025.1',
                        'display': 'Salmonella enterica subsp. enterica serovar Typhimurium str. D23580',
                        'locus_tag_prefix': 'STMMW_'},
}
with open('ESSENTIAL_GENES.json', 'r') as f:
    ESSENTIAL_GENES = json.load(f)

# Study citations for essential gene datasets
ESSENTIAL_GENE_STUDIES = {key: 'Ghomi, Jung et al. 2024' for key in ESSENTIAL_GENES}
ESSENTIAL_GENE_STUDIES['e_coli'] = 'EcoGene'


def _preset_choices():
    """Build the preset choices list used by all forms."""
    return ([('upload', 'Own files'), ('ncbi', 'NCBI assembly accession')]
            + [(key, PRESETS[key]['display']) for key in PRESETS])


def _rewrite_gff_locus_tags(gff_path, prefix):
    """Rewrite GFF locus_tag attributes using old_locus_tag values matching prefix.

    NCBI RefSeq GFFs re-annotate locus tags with *_RS* format. The original tags
    are preserved in old_locus_tag attributes on gene lines (sometimes URL-encoded
    with %2C-separated multiples). CDS/RNA lines typically lack old_locus_tag.

    Two-pass approach:
    1. Build a mapping {new_RS_tag -> old_tag} from lines that have old_locus_tag
    2. Apply the mapping to ALL lines that have a locus_tag
    """
    from urllib.parse import unquote
    with open(gff_path, 'r') as f:
        lines = f.readlines()

    # Pass 1: build mapping from new RS tags to old tags
    tag_map = {}  # e.g. {'T_RS17720': 't3488'}
    for line in lines:
        if line.startswith('#') or '\t' not in line:
            continue
        m_old = re.search(r'old_locus_tag=([^;\n]+)', line)
        if not m_old:
            continue
        m_new = re.search(r'\blocus_tag=([^;\n]+)', line)
        if not m_new:
            continue
        new_tag = m_new.group(1).split(';')[0]  # stop at semicolon if present
        old_tags = unquote(m_old.group(1)).split(',')
        match = next((t for t in old_tags if t.startswith(prefix)), None)
        if match:
            tag_map[new_tag] = match

    # Pass 2: rewrite locus_tag on all feature lines using the mapping
    out = []
    for line in lines:
        if line.startswith('#') or '\t' not in line:
            out.append(line)
            continue
        m_tag = re.search(r'\blocus_tag=([^;\n]+)', line)
        if m_tag:
            current_tag = m_tag.group(1).split(';')[0]
            if current_tag in tag_map:
                line = re.sub(r'\blocus_tag=[^;\n]+', f'locus_tag={tag_map[current_tag]}', line)
        out.append(line)

    with open(gff_path, 'w') as f:
        f.writelines(out)


def _download_ncbi_accession(accession, dest_dir, locus_tag_prefix=None):
    """Download genome FASTA and GFF3 from NCBI using the datasets CLI.

    Returns dict with 'genome' and 'gff' paths.
    If locus_tag_prefix is set, rewrites GFF locus_tag using old_locus_tag values.
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

    if locus_tag_prefix and gff_files:
        _rewrite_gff_locus_tags(gff_files[0], locus_tag_prefix)

    return {'genome': fna_files[0], 'gff': gff_files[0]}


def _convert_genbank_to_fasta_gff(gb_path, dest_dir):
    """Convert a GenBank file to FASTA and GFF3 files.

    Returns dict with 'genome' and 'gff' paths.
    """
    from Bio import SeqIO

    fasta_path = os.path.join(dest_dir, 'genome.fasta')
    gff_path = os.path.join(dest_dir, 'annotation.gff3')

    try:
        records = list(SeqIO.parse(gb_path, 'genbank'))
    except Exception as e:
        raise RuntimeError(
            'Could not parse the uploaded file as GenBank format. '
            'Please make sure it is a valid GenBank/GBFF file. '
            f'Error: {e}')
    if not records:
        raise RuntimeError(
            'No sequence records found in the GenBank file. '
            'Please check that the file is not empty and is in GenBank format (.gb, .gbk, .gbff).')

    # Write FASTA
    SeqIO.write(records, fasta_path, 'fasta')

    # Write GFF3
    with open(gff_path, 'w') as gff:
        gff.write('##gff-version 3\n')
        for record in records:
            seq_id = record.id
            for feature in record.features:
                if feature.type in ('source', 'misc_feature'):
                    continue
                start = int(feature.location.start) + 1  # GFF is 1-based
                end = int(feature.location.end)
                strand = '+' if feature.location.strand == 1 else '-'

                attrs = []
                qualifiers = feature.qualifiers
                if 'locus_tag' in qualifiers:
                    attrs.append(f"locus_tag={qualifiers['locus_tag'][0]}")
                if 'gene' in qualifiers:
                    attrs.append(f"gene={qualifiers['gene'][0]}")
                if 'product' in qualifiers:
                    attrs.append(f"product={qualifiers['product'][0]}")
                if 'old_locus_tag' in qualifiers:
                    attrs.append(f"old_locus_tag={qualifiers['old_locus_tag'][0]}")

                feat_id = qualifiers.get('locus_tag', qualifiers.get('gene', ['unknown']))[0]
                attrs.insert(0, f"ID={feat_id}")

                attr_str = ';'.join(attrs) if attrs else '.'
                gff.write(f"{seq_id}\tGenBank\t{feature.type}\t{start}\t{end}\t.\t{strand}\t.\t{attr_str}\n")

    return {'genome': fasta_path, 'gff': gff_path}


def _resolve_genome_paths(form, time_string):
    """Handle preset-or-upload file logic. Returns dict with 'genome' and 'gff' paths."""
    paths = {}
    if form.presets.data == 'upload':
        base_dir = os.path.join(app.root_path, 'static/data', time_string)
        os.makedirs(base_dir, exist_ok=True)

        if form.upload_format.data == 'genbank':
            gb_file = os.path.join(base_dir, form.genbank.data.filename)
            form.genbank.data.save(gb_file)
            paths = _convert_genbank_to_fasta_gff(gb_file, base_dir)
        else:
            # FASTA + GFF upload
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
        preset = PRESETS[form.presets.data]
        if 'accession' in preset:
            base_dir = os.path.join(app.root_path, 'static/data', time_string)
            paths = _download_ncbi_accession(preset['accession'], base_dir,
                                             preset.get('locus_tag_prefix'))
        else:
            paths['genome'] = preset['genome']
            paths['gff'] = preset['gff']

    # Validate files (skip bundled presets which are known-good)
    is_bundled_preset = (form.presets.data in PRESETS
                         and 'accession' not in PRESETS[form.presets.data])
    if not is_bundled_preset:
        _validate_genome_files(paths)

    return paths


def _validate_genome_files(paths):
    """Validate uploaded FASTA and GFF files before running the pipeline.

    Raises RuntimeError with a user-friendly message if validation fails.
    """
    genome_path = paths['genome']
    gff_path = paths['gff']

    # --- FASTA validation ---
    with open(genome_path) as f:
        first_line = f.readline().strip()
    if not first_line.startswith('>'):
        raise RuntimeError(
            'The uploaded FASTA file does not appear to be in FASTA format. '
            'FASTA files must start with a header line beginning with ">".')

    # Check that FASTA has actual sequence content
    has_sequence = False
    with open(genome_path) as f:
        for line in f:
            if not line.startswith('>') and line.strip():
                has_sequence = True
                break
    if not has_sequence:
        raise RuntimeError(
            'The uploaded FASTA file contains no sequence data. '
            'Please check that the file is not empty or contains only headers.')

    # Collect FASTA sequence IDs
    fasta_ids = set()
    with open(genome_path) as f:
        for line in f:
            if line.startswith('>'):
                fasta_ids.add(line[1:].split()[0])

    # --- GFF validation ---
    with open(gff_path) as f:
        first_line = f.readline().strip()
    # Skip comment lines to find the first feature line
    has_features = False
    has_gene_id = False
    gff_seq_ids = set()
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                raise RuntimeError(
                    'The uploaded GFF file does not appear to be in GFF format. '
                    'GFF files must be tab-separated with 9 columns. '
                    'Make sure you are not uploading a CSV or other format.')
            has_features = True
            gff_seq_ids.add(cols[0])
            if 'locus_tag=' in cols[8] or 'gene=' in cols[8] or 'old_locus_tag=' in cols[8]:
                has_gene_id = True
            # Only need to check a few lines
            if len(gff_seq_ids) > 50:
                break

    if not has_features:
        raise RuntimeError(
            'The uploaded GFF file contains no feature entries. '
            'Please check that the file is a valid GFF/GFF3 annotation file.')

    if not has_gene_id:
        raise RuntimeError(
            'The uploaded GFF file does not contain any gene identifiers. '
            'MASON requires locus_tag, gene, or old_locus_tag attributes in the GFF '
            'to identify genes. Please check the 9th column of your GFF file.')

    # Check that GFF sequence IDs match FASTA headers
    if fasta_ids and gff_seq_ids:
        matching = fasta_ids & gff_seq_ids
        if not matching:
            raise RuntimeError(
                f'The sequence IDs in the GFF file ({", ".join(sorted(gff_seq_ids)[:3])}) '
                f'do not match any sequence IDs in the FASTA file ({", ".join(sorted(fasta_ids)[:3])}). '
                'Please make sure the FASTA and GFF files are from the same genome.')


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
        if "." not in dirs and dirs != "ncbi_dataset":
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
            return render_template("start.html", title="Start", essential=ESSENTIAL_GENES, essential_studies=ESSENTIAL_GENE_STUDIES, form=form, autofill=False)
        base_dir = _init_run_directory(time_string, paths)

        if additional_screen == "essential_genes":
            raw_input = request.form.get('essential_genes_input', '')
            tags = [t.strip() for t in re.split(r'[,\n]+', raw_input) if t.strip()]
            with open(os.path.join(base_dir, "essential_genes.txt"), "w") as ef:
                ef.write("\n".join(tags) + "\n")

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

        _usage_logger.info('%s | mason | %s | %s | %s',
                          datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),
                          request.headers.get('X-Forwarded-For', request.remote_addr).split(',')[0].strip(),
                          form.presets.data, '; '.join(target_genes))

        return redirect(url_for('result', result_id=time_string))
    autofill = request.method == 'GET' and request.path == '/startauto'
    return render_template("start.html", title="Start", essential=ESSENTIAL_GENES, essential_studies=ESSENTIAL_GENE_STUDIES, form=form, autofill=autofill)


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
            return render_template("scrambler.html", title="Scrambler", form=form, essential=ESSENTIAL_GENES, essential_studies=ESSENTIAL_GENE_STUDIES)
        base_dir = _init_run_directory(time_string, paths)

        if additional_screen == "essential_genes":
            raw_input = request.form.get('essential_genes_input', '')
            tags = [t.strip() for t in re.split(r'[,\n]+', raw_input) if t.strip()]
            with open(os.path.join(base_dir, "essential_genes.txt"), "w") as ef:
                ef.write("\n".join(tags) + "\n")

        result_custom_id = form.custom_id.data
        pna_seq = form.seq_input.data.upper()
        pna_file_string = os.path.join(base_dir, "pna_input.fasta")

        with open(pna_file_string, "w") as pna_file:
            pna_file.write(">PNA\n" + pna_seq)

        open(os.path.join(path_parent, "logfile_masonscript.log"), "w").close()

        scrambler_mode = request.form.get('scrambler_mode', 'scramble')
        num_mismatches = form.num_mismatches.data if scrambler_mode == 'mismatch' else '0'

        threading.Thread(target=start_scrambler, name="scramblers",
                         args=[str(path_parent) + "/scrambler.sh", paths['genome'],
                               paths['gff'], time_string, time_string, pna_seq,
                               additional_screen, scrambler_mode, num_mismatches]).start()

        with open(os.path.join(base_dir, "inputs.txt"), "a") as inputfile:
            inputfile.write("\n" + result_custom_id + "\n" + pna_seq)
            inputfile.write("\n" + additional_screen)
            inputfile.write("\n" + scrambler_mode)

        _usage_logger.info('%s | scrambler | %s | %s | %s',
                          datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),
                          request.headers.get('X-Forwarded-For', request.remote_addr).split(',')[0].strip(),
                          form.presets.data, pna_seq)

        return redirect(url_for('result_scrambler', result_id=time_string))
    return render_template("scrambler.html", title="Scrambler", form=form, essential=ESSENTIAL_GENES, essential_studies=ESSENTIAL_GENE_STUDIES)


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
            return render_template("checker.html", title="ASO-Checker", essential=ESSENTIAL_GENES, essential_studies=ESSENTIAL_GENE_STUDIES, form=form)
        base_dir = _init_run_directory(time_string, paths)

        if additional_screen == "essential_genes":
            raw_input = request.form.get('essential_genes_input', '')
            tags = [t.strip() for t in re.split(r'[,\n]+', raw_input) if t.strip()]
            with open(os.path.join(base_dir, "essential_genes.txt"), "w") as ef:
                ef.write("\n".join(tags) + "\n")

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

        _usage_logger.info('%s | checker | %s | %s | %d sequences',
                          datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'),
                          request.headers.get('X-Forwarded-For', request.remote_addr).split(',')[0].strip(),
                          form.presets.data, pna_seq.count('>'))

        return redirect(url_for('result_checker', result_id=time_string))
    return render_template("checker.html", title="ASO-Checker", essential=ESSENTIAL_GENES, essential_studies=ESSENTIAL_GENE_STUDIES, form=form)


@app.route("/result/<result_id>")
def result(result_id):
    ctx = _read_result_context(result_id)
    tgenes = ctx['lines'][2]
    add_screen = ctx['lines'][3].strip()
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
    scrambler_mode = ctx['lines'][4].strip() if len(ctx['lines']) > 4 else "scramble"

    return render_template("scrambler_result.html", title="Result", result=result_id, dir_out=ctx['dir_out'],
                           all_output_dirs=ctx['all_output_dirs'],
                           rfin=ctx['rfinished'], genome_file=ctx['ffile'], gff_file=ctx['gfffile'],
                           custom_id=ctx['custom_id'], pnaseq=pnaseq,
                           time=ctx['time'], errs=ctx['errs'], add_screen=add_screen,
                           scrambler_mode=scrambler_mode)


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

    screen_values = ("none", "human", "microbiome", "essential_genes")
    if last_line in ("yes", "no") and second_last not in screen_values:
        # New format with target_gene and use_ml appended
        use_ml = last_line
        target_gene = second_last
        add_screen = lines[-3].strip() if len(lines) >= 3 else "none"
        pnaseq = "".join(lines[2:-3])
    elif last_line in screen_values:
        add_screen = last_line
        pnaseq = "".join(lines[2:-1])
    else:
        add_screen = "none"
        pnaseq = "".join(lines[2:])

    # Read varna_positions.tsv for VARNA display when target gene is set
    varna_asos = []
    warnings = []
    sd_found = False
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
        sd_path = os.path.join(base_dir, "outputs", "sd_position.txt")
        sd_found = os.path.isfile(sd_path) and os.path.getsize(sd_path) > 0

    return render_template("checker_result.html", title="Result ASO-Checker", result=result_id,
                           dir_out=ctx['dir_out'],
                           all_output_dirs=ctx['all_output_dirs'],
                           rfin=ctx['rfinished'], genome_file=ctx['ffile'], gff_file=ctx['gfffile'],
                           custom_id=ctx['custom_id'], pnaseq=pnaseq,
                           time=ctx['time'], errs=ctx['errs'], add_screen=add_screen,
                           target_gene=target_gene, use_ml=use_ml,
                           varna_asos=varna_asos, warnings=warnings, sd_found=sd_found)


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
