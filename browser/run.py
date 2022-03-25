"""
Copyright (c) 2020 Patrick Pfau
This project is licensed under the terms of the MIT license.

This will start the server
"""

from pnag import app 	# import from the __init__.py in pnag

if __name__ == "__main__":
	app.run(debug=True)  # change after debugging!
