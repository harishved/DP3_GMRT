# DP3_GMRT
Some scripts to process uGMRT data with LOFAR software tools (DP3)

Here I assume that the usual casa based reductions have been applied 
the end product of which is a MS file with only the target scans. 

DP3 and the LOFAR software stack is likely to be useful for direction dependent calibration
and peeling of bright sources, especially in band3 which has a rather large FoV.

So there is no huge advantage in using DP3 for the basic casa based calibration steps with the usual 
primary and secondary calibrator set up. 

The basic steps are to first split the visibility into its different scans
This is because some steps in DP3 have problems dealing with data on a non-uniform time grid

I also noticed that the uGMRT data has small millisecond level offsets from a regular time grid even within a scan.
Turns out that even this is a problem for DP3 so one needs to "force" the data to be on a regular time grid. 
This can be done by brute force rewriting of the TIME column as the ms-level offsets should not have a noticeable error for most imaging applications. 

Then there is the problem of DP3 not liking the fact that my data only has two correlations (RR and LL). It wanted to have full Stokes so I had to change the dimensions of the many columns and write zeros when cross-hand data was unavailable. 

With these changes (all implemented with simple python scripts that use python-casacore), the data can be prepared for DP3.


