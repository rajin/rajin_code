rajin_code
==========
This piece of code is for the spectroscopic reduction of long-slit data from SALT (Southern African Large Telescope).
It is by no means the only way to reduce the data. What it provides though is convenience. Instead of running a few 
IRAF task by hand and on single or list of images, this piece of code, does the reduction at one go.

====================================================================================================================

The Requirements:
This code has been tested to run on IRAF 2.14, Pyraf 1.10 and Python 2.6
There is a know issue with Python 2.7 that IRAF pop-up windows fails to load due to some issue with Tkinter. In some 
cases though it could run.
It would be good to load the following packages: noao, onedspec, twodspec, longslit, apextract, stsdas, in your IRAF 
default login.cl file found in your IRAf home folder.
Also there are some files (ARC line-list, LA-Cosmic file) that you need to paste in your IRAF home folder for the code
to work

====================================================================================================================

Input:
How to input data so that the code run smoothly?
The code can only calibrate 1 ARC spectra when running.
Therefore it would be good idea to select all images (science, flats, etc) to be calibrated by that ARC and paste 
all of them (along with the ARC) in a folder.
Also put the code inside the directory where the images are to run it.
