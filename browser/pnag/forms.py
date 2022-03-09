"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This sets the Forms used on the Website. These variables can be called in the html for display and access.
"""
import os
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import StringField, PasswordField, SubmitField, BooleanField, SelectField, SelectMultipleField
from wtforms.fields import IntegerField
from wtforms.validators import Length, Email, EqualTo, ValidationError, NumberRange, InputRequired, Regexp
from flask_login import current_user
from pnag import app
from pnag.models import User
import re


class NoValidationSelectMultipleField(SelectMultipleField):
	"""
	Overrides the validation of the SelectMultipleField.
	Is necessary because the choices will be added dynamically by JavaScript.
	"""
	def pre_validate(self, form):
		pass


class RegistrationForm(FlaskForm):
	username = StringField('Username', validators=[InputRequired(), Length(min=5, max=20)])
	email = StringField('Email', validators=[InputRequired(), Email()])
	password = PasswordField('Password', validators=[InputRequired(), Length(min=6)])
	confirm_password = PasswordField('Confirm Password', validators=[InputRequired(), EqualTo('password')])
	submit = SubmitField('Sign Up')

	def validate_username(self, username):
		user = User.query.filter_by(username=username.data).first()
		if user:
			raise ValidationError('This username already exist. Please choose a different one.')

	def validate_email(self, email):
		user = User.query.filter_by(email=email.data).first()
		if user:
			raise ValidationError('This email already exist. Please choose a different one.')


class LoginForm(FlaskForm):
	email = StringField('Email', validators=[InputRequired(), Email()])
	password = PasswordField('Password', validators=[InputRequired()])
	remember = BooleanField('Remember me')
	submit = SubmitField('Login')


class UpdateAccountForm(FlaskForm):
	username = StringField('Username', validators=[InputRequired(), Length(min=5, max=20)])
	email = StringField('Email', validators=[InputRequired(), Email()])
	submit = SubmitField('Update')

	def validate_username(self, username):
		if username.data != current_user.username:
			user = User.query.filter_by(username=username.data).first()
			if user:
				raise ValidationError('This username already exist. Please choose a different one.')

	def validate_email(self, email):
		if email.data != current_user.email:
			user = User.query.filter_by(email=email.data).first()
			if user:
				raise ValidationError('This email already exist. Please choose a different one.')


class startForm(FlaskForm):
	custom_id = StringField('Custom ID for result recognition', validators=[InputRequired(), Length(min=3, max=25)])
	genome = FileField('Whole genome FASTA File', validators=[FileAllowed(['fasta', 'fa', 'fna'])])
	gff = FileField('GFF File', validators=[FileAllowed(['gff', 'gff3', 'gff2', 'gff1', 'gtf'])])
	presets = SelectField('Select one of the preset genomes or use "Own files" to upload your own genome', choices=[], validators=[InputRequired()])
	essential = NoValidationSelectMultipleField('Essential genes to be selected', choices=[])
	genes = StringField('Enter the locus tags of target genes separated by comma')
	len_PNA = IntegerField('Length of ASOs', validators=[InputRequired(), NumberRange(min=5, max=19)])
	mismatches = IntegerField('Allowed mismatches for off targets', validators=[InputRequired(), NumberRange(min=0, max=5)])
	bases_before = StringField('Bases before (5 prime) CDS (start codon) to start ASO design (optional)',
							   validators=[Regexp("^[1-9]$|^1[0-9]$|^2[0-5]$|^$", message="Numbers between 1-25 are acepted")])
	submit = SubmitField('Start')

	def validate_genes(self, genes):
		if not genes.data and not self.essential.data:
			raise ValidationError("Please specify any target gene.")
		if genes.data.find('>') != -1:
			raise ValidationError("Don't use the fasta record of the gene. Use the locus_tag instead.")
		elif genes.data.find(';') != -1:
			raise ValidationError("Please separate the locus_tags by comma: ','")

	def validate_presets(self, presets):
		if presets.data == 'upload':
			if not self.genome.data or not self.gff.data:
				raise ValidationError('Please select a preset or upload both files, genome and gff.')


class RequestResetForm(FlaskForm):
	email = StringField('Email', validators=[InputRequired(), Email()])
	submit = SubmitField('Request password reset')

	def validate_email(self, email):
		user = User.query.filter_by(email=email.data).first()
		if user is None:
			raise ValidationError('There is no account with that email. You must register first.')


class ResetPasswordForm(FlaskForm):
	password = PasswordField('Password', validators=[InputRequired()])
	confirm_password = PasswordField('Confirm Password', validators=[InputRequired(), EqualTo('password')])
	submit = SubmitField('Save new password')

