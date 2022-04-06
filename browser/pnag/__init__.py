"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.
"""
from flask import Flask
# for encrypt the password
from flask_bcrypt import Bcrypt
# to sent emails
from flask_mail import Mail
import json

with open('./config.json', 'r') as config_file:
    config = json.load(config_file)

app = Flask(__name__)
# the Key is used by the app to protect from various attacks
app.config['SECRET_KEY'] = config.get('SECRET_KEY')


bcrypt = Bcrypt(app)
# to sent the reset password emails
app.config['MAIL_SERVER'] = "smtp.gmail.com"  #'smtp.web.de'
app.config['MAIL_PORT'] = 587  # 587
app.config['MAIL_USE_TLS'] = True
app.config['MAIL_USE_SSL'] = False
# both need to be set before installation
app.config['MAIL_USERNAME'] = config.get('EMAIL_USER')
app.config['MAIL_PASSWORD'] = config.get('EMAIL_PASS')
mail = Mail(app)

from pnag import routs

