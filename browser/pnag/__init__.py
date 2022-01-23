"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.
"""
from flask import Flask
# for sql database
from flask_sqlalchemy import SQLAlchemy
# for encrypt the password
from flask_bcrypt import Bcrypt
# to handel Login
from flask_login import LoginManager
# to sent emails
from flask_mail import Mail
import json

with open('./config.json', 'r') as config_file:
    config = json.load(config_file)

app = Flask(__name__)
# the Key is used by the app to protect from various attacks
app.config['SECRET_KEY'] = config.get('SECRET_KEY')
# tells the path to the database
app.config['SQLALCHEMY_DATABASE_URI'] = config.get('SQLALCHEMY_DATABASE_URI')
db = SQLAlchemy(app)
bcrypt = Bcrypt(app)
login_manager = LoginManager(app)
login_manager.login_view = 'login'
login_manager.login_message_category = 'info'
# to sent the reset password emails
app.config['MAIL_SERVER'] = 'mail.gmx.net'
app.config['MAIL_PORT'] = 587
app.config['MAIL_USE_TLS'] = True
# both need to be set before installation
app.config['MAIL_USERNAME'] = config.get('EMAIL_USER')
app.config['MAIL_PASSWORD'] = config.get('EMAIL_PASS')
mail = Mail(app)

from pnag import routs

