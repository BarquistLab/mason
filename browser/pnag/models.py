"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This defines the structure of the database with one table for user and one for result.
"""
from pnag import db, login_manager, app
# used for timed tokens
from itsdangerous import TimedJSONWebSignatureSerializer as Serializer
from datetime import datetime
from flask_login import UserMixin
import pickle


@login_manager.user_loader
def load_user(user_id):
	return User.query.get(int(user_id))


class User(db.Model, UserMixin):
	id = db.Column(db.Integer, primary_key=True)
	username = db.Column(db.String(20), unique=True, nullable=False)
	email = db.Column(db.String(120), unique=True, nullable=False)
	password = db.Column(db.String(60), nullable=False)
	# this is a backref to call all results created by this user
	results = db.relationship('Result', backref='owner', lazy=True)

	def get_reset_token(self, expires_sec=600):
		s = Serializer(app.config['SECRET_KEY'], expires_sec)
		return s.dumps({'user_id': self.id}).decode('utf-8')

	@staticmethod
	def verify_reset_token(token):
		s = Serializer(app.config['SECRET_KEY'])
		try:
			user_id = s.loads(token)['user_id']
		except:
			return None
		return User.query.get(user_id)


	def __repr__(self):
		return f"User('{self.username}', '{self.email}')"


class Result(db.Model):
	id = db.Column(db.Integer, primary_key=True)
	custom_id = db.Column(db.String(30), nullable=False)
	genome = db.Column(db.String(35), nullable=False)
	gff = db.Column(db.String(35), nullable=False)
	genes = db.Column(db.String(35), nullable=False)
	mismatches = db.Column(db.Integer, nullable=False)
	finish = db.Column(db.Boolean, nullable=False, default=False)
	result = db.Column(db.String(35))
	pnas = db.Column(db.PickleType)
	date = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
	# this is a backref to call the user who created/owns this result
	user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)

	def __repr__(self):
		return f"Result('{self.result}', '{self.date}')"
