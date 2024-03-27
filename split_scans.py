# Split each unique scan from an MS into a new table
# This uses TAQL
# Intention here is to have tables with TIME on a regular grid that is assumed by DP3
# Note: There is a deep copy here (remove "as plain" is you just want a reference table)
# Author: Harish Vedantham
# 26 March 2024
#
import casacore.tables as tab
import sys
import os
#
ms=sys.argv[1]
t=tab.table(ms)
unique_scans = list(set(list(t.getcol("SCAN_NUMBER"))))
print ("Unique scan numbers = "+str(unique_scans))
nscan = len(unique_scans)
for i in range(nscan):
    newms = ms[:-3]+"_%02d.ms"%i
    os.system("taql select from %s where SCAN_NUMBER=%d giving %s as plain"%(ms,unique_scans[i],newms))


