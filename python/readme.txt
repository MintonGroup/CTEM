1)	Install a distribution of Python 2.x to your machine. This source code
	has been tested on Anaconda Python 4.1.1.

2)  Insure required modules are installed:
	numpy
	os
	scipy
	shutil
	subprocess
	
3)	Place the .py files in this directory into the directory you wish
	to execute from, with the required initialization files.  For this,
	the IDL directory may be copied, and ctem.in replaced by the one found
	in this directory.
	
4)  Inside Anaconda Python, the program can be run through the Spyder2 GUI.
	Outside of Anaconda, the progam can be begun by configuring the ctem.in
	file, then typing
	
	python < ctem_driver.py
	
	at the command line.
