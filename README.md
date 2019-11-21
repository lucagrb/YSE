# YSE
codes for the YSE survey (and possibly more...)

### PS1_query.py

This simple code makes a query to the PS1 catalog and image server and provides three distinct files as output
1) a .fits file (named for the moment as test_run.fits) of the region of the sky, given ra, dec and a width for the bow image
2) a .jpg file (test_run.jpg) showing the same image of the sky with overplotted the sources that are in the PS1 catalog and marked with a blue circle
3) a .csv table (test_run.csv) with the list of all the sources and magnitudes in grizy bands

Of course, this is still a first test.
It is possible to implement the datatab output in other codes, like EAZY, to infer the photo-z

# How to launch

It is simple :) 
From your bash shell type:

> python PS1_query.py ra dec width

Where ra is the Right Ascension of the center of the desired field of view in degrees, dec the Declination of the center of the desired field of view in degrees and width the width of the image (assumed a squared one) in arcmin.

Example - for querying M57, and to obtain the same test_run files that are in this archive, type

> python PS1_query.py 283.395889 33.028572 4

NB i am working to remove some of the warnings that appear during the execution of the code

# PREREQUISITES

Probably you need to install the following python packages

- astroquery
- astropy
- pandas
- photutils

 in addition to the classical python packages like matplotlib and numpy
 
 # CREDITS
 
 Part of this code has been developed starting from other public codes taken from the Astropy community and the PanStarrs collaboration.
 Special thanks go to Renato Martone for the support in providing part of these scripts.
