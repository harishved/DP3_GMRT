# Some utility functions I write when figuring out how to process uGMRT band-2 data
# with DP3. See also https://github.com/harishved/DP3_GMRT/blob/main/gmrt2dp3.py
# Lot of these functions are ad-hoc and were used to troubleshoot
# as I was developing the pipeline
#
import casacore.tables as tab
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import copy
#
#
# Global variable (# of baselines)
# This should be found out from the ms itself but
# for now it is hardcoded
# uGMRT data usually does not have autocorrelations
# and usually 30 antennas are present giving nabse = 435
nbase=435
#
def unflag_all(ms):
    # Unset all flags in a ms
    # can also be done using taql
    # taql update <msname> set FLAG=False
    #
    t = tab.table(ms,readonly=False)
    f = t.getcol("FLAG")
    f[:,:,:]=False
    t.putcol("FLAG",f)
    t.close()

def unflag_cross(ms):
    # Unflag all cross pol products 
    # e.g. RL LR 
    # ASsumes the cross products are in position 1 and 2 (0 indexing)
    t = tab.table(ms,readonly=False)
    f = t.getcol("FLAG")
    f[:,:,1]=False
    f[:,:,2]=False
    t.putcol("FLAG",f)
    t.close()


def getbase(ms,a1,a2):
    # Get the visibilities and flags for a given baseline
    t = tab.table(ms)
    a1list = t.getcol("ANTENNA1",0,nbase,1)
    a2list = t.getcol("ANTENNA2",0,nbase,1)
    I = np.where (np.logical_and(a1list==a1,a2list==a2))[0][0]
    d = t.getcol("DATA",I,t.nrows(),nbase)
    f = t.getcol("FLAG",I,t.nrows(),nbase)
    t.close()
    return d,f

def rfifrac(ms):
    # Return the rfi fraction
    t = tab.table(ms)
    f = t.getcol("FLAG")[:,:,0]
    t.close()
    return (np.sum(f)/np.prod(f.shape))


def plotbase(ms,a1,a2):
    # make a waterfall plot of a given baseline
    # amp and phase
    # Save to plots/bline.pdf
    d1,f1 = getbase(ms,a1,a2)

    plt.figure(figsize=(10,10))
    
    d = d1[:,:,0]
    f = f1[:,:,0]
    da = np.absolute(d)

    plt.subplot(221)
    med = np.median(da[~f])
    mad = np.median(np.absolute(da[~f]-med))
    plt.imshow(np.transpose(da),\
            aspect="auto",origin="lower",\
            interpolation="none",vmin=med-15*mad,\
            vmax=med+15*mad,cmap="viridis")
    plt.subplot(222)
    plt.imshow(np.transpose(np.angle(d)),\
            aspect="auto",origin="lower",\
            interpolation="none",vmin=-np.pi,\
            vmax=np.pi,cmap="hsv")

    if d1.shape[2]==4:
        ii=3
    else:
        ii=1
    d = d1[:,:,ii]
    f = f1[:,:,ii]
    da = np.absolute(d)

    plt.subplot(223)
    med = np.median(da[~f])
    mad = np.median(np.absolute(da[~f]-med))
    plt.imshow(np.transpose(da),\
            aspect="auto",origin="lower",\
            interpolation="none",vmin=med-15*mad,\
            vmax=med+15*mad,cmap="viridis")
    plt.subplot(224)
    plt.imshow(np.transpose(np.angle(d)),\
            aspect="auto",origin="lower",\
            interpolation="none",vmin=-np.pi,\
            vmax=np.pi,cmap="hsv")

    plt.tight_layout()
    plt.savefig("plots/bline.pdf")
    plt.close()
    return


def flag_dropout(msname,nant=30,minlevel=0.5,maxlevel=10):
    # Finding dropout visibility as ones with very low amplitudes
    # Should work on uncalibrated data (hopefully!)
    # For each baseline
    # Find the median amplitude across channels
    # And then find timeslots where the median is lower than some value
    # And flag those timeslots
    nbase = int(nant*(nant-1)/2)
    t = tab.table(msname,readonly=False)
    ntime = int(t.nrows()/nant)
    plt.figure(figsize=(10,6))

    print ("")

    for i in range(nbase):
        d = np.absolute(t.getcol("DATA",i,ntime,nbase))
        f = t.getcol("FLAG",i,ntime,nbase)

        level = np.ones((d.shape[0],))
        for j in range(len(level)):
            level[j] = np.median(d[j,~f[j,:,0],0])

        level/=np.nanmedian(level)
        I = np.where(np.logical_or(level>maxlevel,level<minlevel))[0]

        for ii in I:
            f[ii,:,0]=True
        plt.plot(level,'g',linewidth=0.1)

        level = np.ones((d.shape[0],))
        for j in range(len(level)):
            level[j] = np.median(d[j,~f[j,:,1],1])

        level/=np.nanmedian(level)
        I = np.where(np.logical_or(level>maxlevel,level<minlevel))[0]
   
        for ii in I:
            f[ii,:,1]=True
            
        #t.putcol("FLAG",f,i,ntime,nbase)
        plt.plot(level,'m',linewidth=0.1) 
        print ("\r %d/%d"%(i,nbase))

    t.close()
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig("plots/eraseme.pdf")
    plt.close()

def flag_filter(ms):
    # ASsumes 8192 channels
    # Make it form for other nchan in the future
    t = tab.table(ms,readonly=False)
    f = t.getcol("FLAG")
    f[:,0:2300,:]=True
    f[:,4500:5600,:]=True
    f[:,7000:,:]=True
    t.putcol("FLAG",f)
    t.close()


def smooth_cal_phase(d):
   # Remove RFI prone channles in phase solution
   # Then find a single phase solution over time
   # as a 2nd order polynomial fit to remaining phase
   # This will fail if there are phase wraps
   # But this simple polyfit seems to work as the residuals delays appear small
   #
   ntime,nfreq,nant,npol = d.shape
   part = int(np.round(nfreq/294*125))

   f = np.arange(0,nfreq,1.0)
   stat = np.zeros((nfreq,))

   for i in range(nant):
      for j in range(npol):
         stat += np.absolute(np.nanmean(np.exp(1j*d[:,:,i,j]),axis=0))
         I = np.where(stat<0.75)

   I = np.where(stat<0.97)[0]
   d[:,I,:,:]=np.nan

   dnew = np.zeros(d.shape)

   plt.figure(figsize=(12,6))
   for i in range(nant):
      for j in range(npol):
         stat = np.nanmean(d[:,:,i,j],axis=0)
         plt.plot(stat,'k',linewidth=0.2)

         x = f[:part]
         y = stat[:part]
         I = ~np.isnan(y)

         p = np.polyfit(x[I],y[I],2)
         pv = np.polyval(p,f[:part])
         for k in range(ntime):
            dnew[k,:part,i,j] = pv

         plt.plot(f[:part],pv,'m',linewidth=0.4)

         x = f[part:]
         y = stat[part:]
         I = ~np.isnan(y)

         p = np.polyfit(x[I],y[I],2)
         pv = np.polyval(p,f[part:])
         for k in range(ntime):
            dnew[k,part:,i,j] = pv

         plt.plot(f[part:],pv,'m',linewidth=0.4)

   plt.xlabel("Freq index")
   plt.ylabel("Phase in radian")
   plt.minorticks_on()
   plt.tight_layout()
   plt.savefig("plots/phase_smooth.pdf")
   plt.close()

   return dnew

def smooth_cal_amp(oldamp):


   ntime,nfreq,nant,npol = oldamp.shape
   part = int(np.round(nfreq/294*125))

   newamp = oldamp.copy()

   badchans = []
   badchans.append([0,17])
   badchans.append([19,21])
   badchans.append([70,100])
   clip_up = 100
   clip_down = 0.1
   print (newamp.shape)
   print ("\n\n\n - - - - - - - - -  HERE - - - - - - - - - \n\n\n")
   newamp[newamp>clip_up]=np.nan

   newamp[newamp<clip_down]=np.nan

   for x in badchans:
      newamp[:,int(x[0]/294.*nfreq):int(x[1]/294.*nfreq)+1,:,:] = np.nan

   freq = np.arange(0.0,nfreq,1.0)

   for i in range(nant):
      for j in range(npol):
         stat = np.nanmedian(newamp[:,:,i,j],axis=0)

         x = freq[:part]
         y = stat[:part]
         I = ~np.isnan(y)

         try:
            p = np.polyfit(x[I],y[I],1)
         except:
            print ("Polyfit failed massively, Taking the median amplitude instead")
            if np.sum(I)==0:
               print ("It is due to an all NaN slice. Assuming unity gain")
               p = [0,0,1]
            else:
               np.savez("polyfit_fail_%s.npz"%(datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")),\
                        iant=i,ipol=j,freq=x[I],amp=y[I])
               p = [0,0,np.nanmedian(y[I])]
         pv = np.polyval(p,freq[:part])
         for k in range(ntime):
            newamp[k,:part,i,j]= pv

         x = freq[part:]
         y = stat[part:]
         I = ~np.isnan(y)
         try:
            p = np.polyfit(x[I],y[I],2)
         except:
            print ("Polyfit failed massively, Taking the median amplitude instead")
            if np.sum(I)==0:
               print ("It is due to all NaN slice. ASsuming unity gain")
               p = [0,0,1]
            else:
               np.savez("polyfit_fail_%s.npz"%(datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")),\
                        iant=i,ipol=j,freq=x[I],amp=y[I])
               p = [0,0,np.nanmedian(y[I])]

         pv = np.polyval(p,freq[part:])
         for k in range(ntime):
            newamp[k,part:,i,j]= pv

   plt.figure(figsize=(12,6))
   for i in range(nant):
      for j in range(npol):
         oldstat = np.nanmedian(oldamp[:,:,i,j],axis=0)
         newstat = np.nanmedian(newamp[:,:,i,j],axis=0)
         plt.plot(oldstat,'k',linewidth=0.2)
         plt.plot(newstat,'m',linewidth=0.2)
   plt.xlabel("Channel index")
   plt.ylabel("Gain amplitude")
   plt.minorticks_on()
   plt.ylim([0.2,10])
   plt.yscale("log") 
   plt.close()


   return newamp


def amp_v_blen(msname):
   # Plot the vis amplitude versus uv distance
   # Included all unflagged data
   #
   t = tab.table(msname)
   d = t.getcol("DATA")
   f = t.getcol("FLAG")
   uvw = t.getcol("UVW")
   blen = np.sqrt(np.sum(uvw**2,axis=1))
   t.close()
   t = tab.table(msname+"/SPECTRAL_WINDOW")
   freq = t.getcol("CHAN_FREQ").flatten()
   t.close()
   lam = 3e8/freq

   plt.figure(figsize=(12,12))
   plt.subplot(211)
   for i in range(len(freq)):
      ff = ~f[:,i,0]
      plt.plot(blen[ff]/lam[i],np.absolute(d[ff,i,0]),'m.',ms=0.02)
   plt.yscale("log")
   plt.minorticks_on()
   plt.subplot(212)
   for i in range(len(freq)):
      ff = ~f[:,i,1]
      plt.plot(blen[ff]/lam[i],np.absolute(d[ff,i,1]),'g.',ms=0.02)
   plt.minorticks_on()
   plt.yscale("log")
   plt.xlabel("uv distance")
   plt.ylabel("Vis Amplitude / Jy")
   plt.tight_layout()
   plt.savefig("plots/ampvdist.png")
   plt.close()


