"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.
"""
from flask import Flask
import json

with open('./config.json', 'r') as config_file:
    config = json.load(config_file)

app = Flask(__name__)
# the Key is used by the app to protect from various attacks
app.config['SECRET_KEY'] = config.get('SECRET_KEY')

from pnag import routs
