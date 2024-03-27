# Force the time column to be on a regular grid
# GMRT data had irregularities at the ~1 millisecond timescale
# And DP3 misbehaves because of this
# So this code just forces the time values to be on a regular grid
# This should only lead to ms level timing errors which should not matter for most
# imaging applications. 
# Author: Harish Vedantham
# 27 March 2024
#
import casacore.tables as tab
import numpy as np
import sys
#
# Check nant and nbase (wrong values here will lead to problems!)
nant=30
nbase=int(nant*(nant-1)/2)
#
ms=sys.argv[1]
t=tab.table(ms,readonly=False) # Open the MS for read/write

old_times = t.getcol("TIME") # Existing time data
nt=int(len(old_times)/nbase) # Number of integrations
#
old_times_unique=np.reshape(old_times,(nt,nbase))[:,0] # List of timestamps
dt=(old_times_unique[-1]-old_times_unique[0])/(nt-1) # Time interval
#
print ("Mean interval = %f"%dt)
new_times_unique = np.arange(nt)*dt+old_times_unique[0] # New list of timestamps
new_times=np.zeros(old_times.shape,dtype=type(old_times[0])) # Initialize new time columns data
# Fill in new time column data
for i in range(nt):
    new_times[i*nbase:i*nbase+nbase]=new_times_unique[i]
#
# Optionally rename old time column and write new data into the TIME column
# Not necessary perhaps as the timing error is very small and should lead to 
# insignificant errors. So its fine to just overvrite the TIME column
'''
coldesc = t.getcoldesc("TIME") # Column descriptor
dminfo = t.getdminfo("TIME")
t.renamecol("TIME","OLD_TIME") # Rename existing time column
t.flush()
t.addcols(tab.maketabdesc(tab.makecoldesc("TIME",coldesc)),dminfo) # Add new time column
'''
t.putcol("TIME",new_times) # copy data into new time column
interval=np.ones(new_times.shape)*dt
t.putcol("INTERVAL",interval)
exposure=np.ones(new_times.shape)*dt
t.putcol("INTERVAL",exposure)
t.putcol("TIME_CENTROID",new_times)
t.flush()
t.close()
