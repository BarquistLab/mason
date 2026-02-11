"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This sets the Forms used on the Website. These variables can be called in the html for display and access.
"""
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import StringField, SubmitField, SelectField, SelectMultipleField, BooleanField
from wtforms.widgets import TextArea
from wtforms.fields import IntegerField
from wtforms.validators import Length, ValidationError, NumberRange, InputRequired, Regexp, Optional
from pnag import app


class NoValidationSelectMultipleField(SelectField):
    """
	Overrides the validation of the SelectMultipleField.
	Is necessary because the choices will be added dynamically by JavaScript.
	"""

    def pre_validate(self, form):
        pass


class BaseOrganismForm(FlaskForm):
    """Base form with fields shared by all tools: custom ID, genome/gff upload, and presets."""
    custom_id = StringField('Custom ID for result recognition', validators=[InputRequired(), Length(min=3, max=25)])
    genome = FileField('FASTA File', validators=[FileAllowed(['fasta', 'fa', 'fna'])])
    gff = FileField('GFF File', validators=[FileAllowed(['gff', 'gff3', 'gff2', 'gff1', 'gtf'])])
    presets = SelectField('Genome of target organism', choices=[], validators=[InputRequired()])
    ncbi_accession = StringField('NCBI Assembly Accession',
                                 validators=[Optional(), Regexp(r'^GC[AF]_\d{9}(\.\d+)?$',
                                             message='Enter a valid NCBI assembly accession (e.g. GCF_000005845.2)')])

    def validate_presets(self, presets):
        if presets.data == 'upload':
            if not self.genome.data or not self.gff.data:
                raise ValidationError('Please select a preset or upload both files, genome and gff.')
        elif presets.data == 'ncbi':
            if not self.ncbi_accession.data or not self.ncbi_accession.data.strip():
                raise ValidationError('Please enter an NCBI assembly accession number.')


class startForm(BaseOrganismForm):
    essential = NoValidationSelectMultipleField('Select essential genes', choices=[])
    genes = StringField('Enter the locus tags of target genes separated by comma',
                        validators=[Regexp(
                            "^[^,;]*$|^[^,;]+,[^,;]+$|^[^,;]+(,[^,;]+){2}$|^[^,;]+(,[^,;]+){3}$|^[^,;]+(,[^,;]+){4}$",
                            message="please use less than 5 genes")])
    len_PNA = IntegerField('Length of ASOs', validators=[InputRequired(), NumberRange(min=7, max=16)])

    bases_before = StringField('Bases before (5 prime) CDS (start codon) to start ASO design (optional)',
                               validators=[
                                   Regexp("^[1-9]$|^1[0-9]$|^2[0]$|^$", message="Numbers between 1-20 are acepted")])
    use_ml = BooleanField('Use machine-learning (random forest) model to predict MICs', default=False)
    submit = SubmitField('Submit & start MASON')

    def validate_genes(self, genes):
        if not genes.data and not self.essential.data:
            raise ValidationError("Please specify any target gene.")
        if genes.data.find('>') != -1:
            raise ValidationError("Don't use the fasta record of the gene. Use the locus_tag instead.")
        elif genes.data.find(';') != -1:
            raise ValidationError("Please separate the locus_tags by comma: ','")


class ScrambledForm(BaseOrganismForm):
    seq_input = StringField("ASO-sequence (5' to 3'):",
                            validators=[InputRequired(), Regexp("^[ATGCatgc]+$", message="use DNA alphabet only (A,T,G,C)"),
                                        Length(min=8, max=15, message="choose sequence with 8-14 nucleotides!")])
    submit_scr = SubmitField('Submit & start Scrambler')


class CheckerForm(BaseOrganismForm):
    seq_input = StringField("ASO-sequence(s) (5' to 3'):", widget=TextArea(),
                            validators=[InputRequired()])
    target_gene = StringField('Target gene locus tag (optional)',
                              validators=[Optional(), Length(max=50)])
    use_ml = BooleanField('Use machine-learning (random forest) model to predict MICs', default=False)
    submit_scr = SubmitField('Submit & start Checker')
