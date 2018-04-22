From:	"Keith Bechtol" <KBechtol@lsst.org>
Date:  	Wed, May 24, 2017 4:51 pm

Hi Rick, cc. Zoheyr,

Yes, the FITS files are far superior! Much faster analysis in 
addition to easier file storage.

I attach three files.

(1) des_io.py - a modified version of the file you sent me from Zoheyr 
that can also parse the SNANA dump files. I needed this for my studies 
because I wanted to know the detection efficiency as a function of 
various event parameters, and it was import to load the full dump file 
at once. Also, this allowed me to select specific events for analysis.

(2) snana_fitsio.py - this is the code that parses the photometry files
in FITS format. Probably the most useful function is selectCID, which 
grabs the photometry information for a single event according to its CID.
To do this, you must pass in the list file. The function searches through
the set of header file to determine which photometry file contains the 
desired event. It then opens the relevant photometry file and returns 
the desired event.

If you want to detect events en masse, there are also functions which do 
this task more efficiently by processing many events in parallel. See the 
detectChuck function, which processes a given photometry file in parallel,
and the detectAll function, which performs a search on all the photometry 
files sequentially.

(3) sensitivity.py - example analysis script that shows how the functions
in snana_fitsio.py might be used. This is what I used for the IceCube 
forecasts.

Please let me know if you have any questions,

Keith

