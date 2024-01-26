"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This sets the Forms used on the Website. These variables can be called in the html for display and access.
"""
import os
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import StringField, PasswordField, SubmitField, BooleanField, SelectField, SelectMultipleField, RadioField, \
    TextAreaField
from wtforms.fields import IntegerField
from wtforms.validators import Length, Email, EqualTo, ValidationError, NumberRange, InputRequired, Regexp
from flask_login import current_user
from pnag import app
import re


class NoValidationSelectMultipleField(SelectField):
    """
	Overrides the validation of the SelectMultipleField.
	Is necessary because the choices will be added dynamically by JavaScript.
	"""

    def pre_validate(self, form):
        pass


class startForm(FlaskForm):
    custom_id = StringField('Custom ID for result recognition', validators=[InputRequired(), Length(min=3, max=25)])
    genome = FileField('FASTA File', validators=[FileAllowed(['fasta', 'fa', 'fna'])])
    gff = FileField('GFF File', validators=[FileAllowed(['gff', 'gff3', 'gff2', 'gff1', 'gtf'])])
    genome_2 = FileField('FASTA File', validators=[FileAllowed(['fasta', 'fa', 'fna'])])
    gff_2 = FileField('GFF File', validators=[FileAllowed(['gff', 'gff3', 'gff2', 'gff1', 'gtf'])])
    presets = SelectField('Genome of target organism', choices=[], validators=[InputRequired()])
    essential = NoValidationSelectMultipleField('Select essential genes', choices=[])
    genes = StringField('Enter the locus tags of target genes separated by comma',
                        validators=[Regexp(
                            "^[^,;]*$|^[^,;]+,[^,;]+$|^[^,;]+(,[^,;]+){2}$|^[^,;]+(,[^,;]+){3}$|^[^,;]+(,[^,;]+){4}$",
                            message="please use less than 5 genes")])
    len_PNA = IntegerField('Length of ASOs', validators=[InputRequired(), NumberRange(min=7, max=16)])
    mismatches = IntegerField('Allowed mismatches for off targets',
                              validators=[InputRequired(), NumberRange(min=0, max=4)])
    bases_before = StringField('Bases before (5 prime) CDS (start codon) to start ASO design (optional)',
                               validators=[
                                   Regexp("^[1-9]$|^1[0-9]$|^2[0]$|^$", message="Numbers between 1-20 are acepted")])
    submit = SubmitField('Submit & start MASON')

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


class ScrambledForm(FlaskForm):
    custom_id = StringField('Custom ID for result recognition', validators=[InputRequired(), Length(min=3, max=25)])
    genome = FileField('FASTA File', validators=[FileAllowed(['fasta', 'fa', 'fna'])])
    gff = FileField('GFF File', validators=[FileAllowed(['gff', 'gff3', 'gff2', 'gff1', 'gtf'])])
    seq_input = StringField("ASO-sequence (5' to 3'):",
                            validators=[InputRequired(), Regexp("^[ATGC]+", message="use DNA alphabet only (A,T,G,C)"),
                                        Length(min=8, max=14, message="choose sequence with 8-14 nucleotides!")])
    presets = SelectField('Genome of target organism', choices=[], validators=[InputRequired()])
    submit_scr = SubmitField('Submit & start Scrambler')
    #submit = SubmitField('Submit & start MASON')

    def validate_presets(self, presets):
        if presets.data == 'upload':
            if not self.genome.data or not self.gff.data:
                raise ValidationError('Please select a preset or upload both files, genome and gff.')



class startautoForm(FlaskForm):
    custom_id = StringField('Custom ID for result recognition', validators=[InputRequired(), Length(min=3, max=25)],
                            default="test")
    genome = FileField('FASTA File', validators=[FileAllowed(['fasta', 'fa', 'fna'])], default="ABC")
    gff = FileField('GFF File', validators=[FileAllowed(['gff', 'gff3', 'gff2', 'gff1', 'gtf'])])
    genome_2 = FileField('FASTA File', validators=[FileAllowed(['fasta', 'fa', 'fna'])])
    gff_2 = FileField('GFF File', validators=[FileAllowed(['gff', 'gff3', 'gff2', 'gff1', 'gtf'])])
    presets = SelectField('Genome of target organism', choices=[], validators=[InputRequired()])
    essential = NoValidationSelectMultipleField('Select essential genes', choices=[])
    genes = StringField('Enter the locus tags of target genes separated by comma',
                        validators=[Regexp(
                            "^[^,;]*$|^[^,;]+,[^,;]+$|^[^,;]+(,[^,;]+){2}$|^[^,;]+(,[^,;]+){3}$|^[^,;]+(,[^,;]+){4}$",
                            message="please use less than 5 genes")], default="b0184")
    len_PNA = IntegerField('Length of ASOs', validators=[InputRequired(), NumberRange(min=7, max=16)], default=10)
    mismatches = IntegerField('Allowed mismatches for off targets',
                              validators=[InputRequired(), NumberRange(min=0, max=4)],
                              default=2)
    bases_before = StringField('Bases before (5 prime) CDS (start codon) to start ASO design (optional)',
                               validators=[
                                   Regexp("^[1-9]$|^1[0-9]$|^2[0]$|^$", message="Numbers between 1-20 are acepted")])
    submit = SubmitField('Submit & start MASON')

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
