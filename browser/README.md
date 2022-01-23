# Webserver
This is a basic Webserver for using the PNA-Generator

## Install
This Webserver needs at least Python3 to run.
### Needed python packages
To install the needed python packages run the following commands in the terminal, when you use linux. You need to be in root folder of this repository.

	sudo apt-get install python3-pip
	pip3 install -r requirements.txt
	
## Start the Server
To start the server just type the following commands in the terminal while being in the directory of the server. 
You need to be in this folder (browser_design).

	python3
	>>>from pnag import db
	>>>db.create_all()
	>>>exit()
	python3 run.py
	

## Preview of result page
![Result Image](Result%20Preview.png)

This project is licensed under the terms of the MIT license.
