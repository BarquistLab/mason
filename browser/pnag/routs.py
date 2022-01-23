"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This defines the structure of the website with its routs and functions within.
It tells which template to use and validates submitted forms.
"""
import os
from datetime import datetime
from flask import render_template, url_for, flash, redirect, request, send_file, Markup
# import needed things from other files in this package
from pnag import app, db, bcrypt, mail
from pnag.forms import (RegistrationForm, LoginForm, UpdateAccountForm, startForm,
                        RequestResetForm, ResetPasswordForm)
from pnag.models import User, Result
from flask_login import login_user, current_user, logout_user, login_required
from flask_mail import Message
import base64
import json
import threading
from start import start_calculation
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
    if current_user.is_authenticated:
        # changes what is shown when the user is logged in
        results = current_user.results
        return render_template("home.html", title="Home", results=reversed(results))
    else:
        return render_template("home.html", title="Home")


@app.route("/register", methods=['GET', 'POST'])
def register():
    if current_user.is_authenticated:
        return redirect(url_for('home'))
    form = RegistrationForm()
    if form.validate_on_submit():
        # secure password to store in db
        hashed_password = bcrypt.generate_password_hash(form.password.data).decode('utf-8')
        email = form.email.data
        user = User(username=form.username.data, email=email.lower(), password=hashed_password)
        db.session.add(user)
        db.session.commit()
        # sets a nice info message
        flash('Your account has been created! You are now able to log in.', 'success')
        return redirect(url_for('login'))
    return render_template("registration.html", title="Register", form=form)


@app.route("/login", methods=['GET', 'POST'])
def login():
    if current_user.is_authenticated:
        return redirect(url_for('home'))
    form = LoginForm()
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        if user and bcrypt.check_password_hash(user.password, form.password.data):
            login_user(user, remember=form.remember.data)
            ''' if the logged in user have tried to access a page with required login and was redirected here, 
			it redirects back to the page which was tried to accessed before'''
            next_page = request.args.get('next')
            return redirect(next_page) if next_page else redirect(url_for('home'))
        else:
            flash('Login unsuccessful. Please check email and password.', 'danger')
    return render_template("login.html", title="Login", form=form)


@app.route("/logout")
def logout():
    logout_user()
    return redirect(url_for('home'))


@app.route("/account", methods=['GET', 'POST'])
@login_required  # site is only accessable when logged in
def account():
    form = UpdateAccountForm()
    if form.validate_on_submit():
        current_user.username = form.username.data
        current_user.email = form.email.data
        db.session.commit()
        flash('Your account has been updated!', 'success')
        return redirect(url_for('account'))
    elif request.method == "GET":
        form.username.data = current_user.username
        form.email.data = current_user.email
    return render_template("account.html", title="Account", form=form)


@app.route("/start", methods=['GET', 'POST'])
@login_required
def start():
    # form to upload the input files
    form = startForm()
    choices = [(key, PRESETS[key][2]) for key in PRESETS.keys()]
    choices.append(('upload', 'Own files'))
    form.presets.choices = choices
    if form.validate_on_submit():
        paths = {}
        time_string = datetime.utcnow().strftime('_%Y_%m_%d_%H_%M_%S_')
        # save uploaded data with timestring attached:
        if form.presets.data == 'upload':
            paths['genome'] = save_file(form.genome.data, 'genome', time_string)
            paths['gff'] = save_file(form.gff.data, 'gff', time_string)
        else:
            # genome first value in list
            files = PRESETS[form.presets.data]
            paths['genome'] = files[0]
            paths['gff'] = files[1]
        target_genes = form.genes.data
        print(target_genes)
        for locus_tag in form.essential.data:
            if len(target_genes) > 1:
                target_genes += ', ' + locus_tag  # add comma only if already a lt was added:
            else:
                target_genes += locus_tag
        paths['genes'] = save_file(target_genes, 'genes', time_string)
        paths['pnas'] = save_file(None, 'pnas', time_string)

        r = Result(custom_id=form.custom_id.data, genome=paths['genome'], genes=paths['genes'], gff=paths['gff'],
                   mismatches=form.mismatches.data, finish=False, result=paths['pnas'], user_id=current_user.id)
        db.session.add(r)
        db.session.commit()
        # Now run MASON as background process while continuing with start.html and showing the "waiting" html:
        threading.Thread(target=start_calculation, name="masons", args=[path_parent.__str__() + "/mason.sh", paths['genome'],
                                                                        paths['gff'], target_genes, str(form.len_PNA.data),
                                                                        str(form.mismatches.data),
                                                                        str(r.id) + "_" + str(r.user_id), r.id]).start()
        return redirect(url_for('result', result_id=r.id))
    return render_template("start.html", title="Start", essential=ESSENTIAL_GENES, form=form)


@app.route("/result/<int:result_id>")
@login_required
def result(result_id):
    # each result gets its own page to access it. Just by calling .../result/<id of result>.
    res = Result.query.get_or_404(result_id)
    dir_out = "../static/data/" + str(res.id) + "_" + str(res.user_id) + "/outputs"
    svg_res1 = dir_out + "/heatmap.png"
    svg_res2 = dir_out + "/plot_ots_whole_transcriptome.png"
    svg_res3 = dir_out + "/tm.png"
    print(svg_res1)
    # just the owner can see the result page of his results
    if current_user == res.owner or current_user.username in ["PatrickPfau", "jakobjung"]:
        return render_template("result.html", title="Result", result=res, svg_res1=Markup(svg_res1),
                               svg_res2=Markup(svg_res2), svg_res3=Markup(svg_res3))
    else:
        return redirect(url_for('home'))


@app.route("/delete_result/<result_id>")
@login_required
def delete_result(result_id):
    res = Result.query.get_or_404(result_id)
    # just the owner can delete the result
    if (current_user == res.owner and res.finish) or current_user.username == "PatrickPfau":
        first = res.genome.split('/')
        if 'presets' not in first:
            os.remove(res.gff)
            os.remove(res.genome)
        os.remove(res.genes)
        try:
            os.remove(res.result)
        except FileNotFoundError:
            pass
        db.session.delete(res)
        db.session.commit()
    else:
        flash('You cannot delete this result.', 'error')
    return redirect(url_for('home'))


@app.route("/download/<path>")
@login_required
def download(path):
    # just for download the file in path
    path = url_decode(path)
    first, f_ext = os.path.splitext(path)
    preset = first
    preset = preset.split('/')
    first.split('_')
    if 'presets' in preset:
        return send_file(path, as_attachment=True)
    elif current_user.id == int(first[-1]):
        # just the owner can download the file
        try:
            return send_file(path, as_attachment=True)
        except FileNotFoundError:
            flash(u'This file does not exist. Maybe no target genes found to build PNA.fasta.', 'danger')
            return redirect(url_for('home'))
    else:
        flash('You need to be the owner of the file.', 'warning')
        return redirect(url_for('home'))


@app.route("/about")
def about():
    return render_template("about.html", title="About")


def send_reset_email(user):
    token = user.get_reset_token()
    msg = Message('Password Reset Request', sender='pna-generator@gmx.de', recipients=[user.email])
    msg.body = f'''To reset your password, visit the following link:
	http://pna-generator.helmholtz-hiri.de{url_for('reset_request_token', token=token)}

If you did not make this request then simply ignore this e-mail and no changes will be made.
This e-mail was created automatically. Please do not reply to this e-mail.
	'''
    mail.send(msg)


@app.route("/reset_password", methods=['GET', 'POST'])
def reset_request():
    if current_user.is_authenticated:
        return redirect(url_for('home'))
    form = RequestResetForm()
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        send_reset_email(user)
        flash('An e-mail has been sent to reset your password. Check your spam folder as well.', 'info')
        return redirect(url_for('login'))
    return render_template('reset_request.html', title='Reset Passowrd', form=form)


@app.route("/reset_password/<token>", methods=['GET', 'POST'])
def reset_request_token(token):
    if current_user.is_authenticated:
        return redirect(url_for('home'))
    user = User.verify_reset_token(token)
    if user is None:
        flash('That is an invalid or expired token', 'warning')
        return redirect(url_for('reset_request'))
    form = ResetPasswordForm()
    if form.validate_on_submit():
        # secure password to store in db
        hashed_password = bcrypt.generate_password_hash(form.password.data).decode('utf-8')
        user.password = hashed_password
        db.session.commit()
        # sets a nice info message
        flash('Your password has been updated! You are now able to log in.', 'success')
        return redirect(url_for('login'))
    return render_template('reset_request_token.html', title='Reset Passowrd', form=form)


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
        file_path = os.path.join(app.root_path, 'temp', file_fn)
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
