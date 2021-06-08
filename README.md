# FCS-Time-Symmetric-Analysis
Contains original Igor Pro code from Ishii &amp; Tahara, Optics Express 2015, of the beautifully elegant time symmetric correlation function that removes after pulsing, together with my Matlab port.

# Correction of the afterpulsing effect in fluorescence correlation spectroscopy using time symmetry analysis
## Kunihiko Ishii and Tahei Tahara
Optics Express Vol. 23, Issue 25, pp. 32387-32400 (2015) •https://doi.org/10.1364/OE.23.032387
https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-23-25-32387&id=333299

## Original Ishii & Tahara Igor Pro code
Igor Pro from https://www.wavemetrics.com/ is a wonderful data analysis package that I really love.  It's a great price and super flexible, but a little esoteric, and the reason I've decide to port this great code from Ishii & Tahara to Matlab and eventually Python

The OSA2015_Ishii_Tahara folder contains the original code that accompanies the paper above.  It runs 32bit Igor Pro, and utilizes a C function via Igor Pro XOP.  Example experiment contains macro time data (in units of laser/sync pulses), micro time (in units of TCSPC bins), decay data (histogram of all macro counts), and cx (tau offset for correlation function; note: 1 length difference to output correlation).  Run `correctedFCS()`; assumes the above waves are in place.

# 
