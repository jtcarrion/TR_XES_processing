import os
import psana
from xtcav2.LasingOnCharacterization import *
import argparse
import h5py
import numpy as np
from scipy import ndimage


parser = argparse.ArgumentParser(description='Batch script for processing Epix100 raw data into VDS structures')
parser.add_argument('-e', dest='exp', type=str, required=True, help='experiment string')
parser.add_argument('-r', dest='run', type=int, required=True, help='run number to process')
parser.add_argument('-o',dest='outdir', type=str, required=True,help='output directory')
parser.add_argument('--laser-on-code',dest='laser_on_code',type=int,default=-1,help='laser on EVR code (ask LCLS staff)')
parser.add_argument('--store-raw',dest='store_raw', action='store_true', help='store the raw images (no pedastal , gain, or geom correction)')
parser.add_argument('--max-events',dest='max_events', type=int, default=np.inf,help='how many events to process... dont use with mpi, -1 for all')
parser.add_argument('-t',dest='tag', type=str, default="epix",help='output tag')
parser.add_argument('--rowmin',type=int, default=0,help='row min for xes integration, default=0')
parser.add_argument('--rowmax',type=int, default=-1,help='row max for xes integration, default = -1')
parser.add_argument('--rotation',type=float, default=-89,help='rotation of image in degrees, default=-89')
parser.add_argument('--threshold',type=float, default=10,help='pixel threshold value, default=10')
parser.add_argument('--eventsperfile',type=int, default=10000,help='events per file, default = 10000')
parser.add_argument('--no-xtcav',dest='no_xtcav',action='store_true',help='ignore xtcav (becuase it is broken...)')
args = parser.parse_args()

# Name/nomenclature of the experiment. (like mfxls0816 for MFX 2016 experiment of Petra)
EXP_STRING=args.exp
subfile=0

def nexth5(args, rank):
    global subfile,h5out,epix,spectra,mtime,mtimenano,fiducial,epixraw, power,pulse,agreement,dsidx,index,laser
    outname = "r%04d-%s-c%03d-r%03d.h5" %(args.run, args.tag,subfile,rank)
    outputpath = os.path.join( args.outdir , outname )
    # small data engine
    h5out=h5py.File(outputpath,'w')
    index = h5out.create_dataset('index', (nf,), dtype='i', chunks=(1,), maxshape=(None,)) 
    laser = h5out.create_dataset('laser',(nf,), dtype='i', chunks=(1,), maxshape=(None,))
    spectra = h5out.create_dataset('spectra',(nf,ns), dtype='f4', chunks=(1,ns), maxshape=(None,ns))
    mtime = h5out.create_dataset('machineTime',(nf,), dtype='f4',chunks=(1,),maxshape=(None,))
    mtimenano = h5out.create_dataset('machineTimeNano',(nf,), dtype='f4',chunks=(1,),maxshape=(None,))
    fiducial = h5out.create_dataset('fiducials',(nf,), dtype='f4',chunks=(1,),maxshape=(None,))
    if not args.no_xtcav:
        power = h5out.create_dataset('xtcav/power',(nf,), dtype='f4',chunks=(1,),maxshape=(None,))
        pulse = h5out.create_dataset('xtcav/pulseFWHM',(nf,), dtype='f4',chunks=(1,),maxshape=(None,))
        agreement = h5out.create_dataset('xtcav/agreement',(nf,), dtype='f4',chunks=(1,),maxshape=(None,))
    if args.store_raw:
        epix = h5out.create_dataset('epix',(nf,n0,n1), dtype='f4', chunks=(1,n0,n1), maxshape=(None,n0,n1))
        epixraw = h5out.create_dataset('epixraw',(nf,n0,n1), dtype='f4', chunks=(1,n0,n1), maxshape=(None,n0,n1))
    dsidx = 0

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#psana
#ds_string = 'exp=%s:run=%d:dir=/cds/data/psdm/cxi/%s/xtc' %(EXP_STRING,  args.run, EXP_STRING) # experiment string
ds_string = 'exp=%s:run=%d:dir=/cds/data/psdm/mfx/%s/xtc' %(EXP_STRING,  args.run, EXP_STRING) # experiment string
print(ds_string)
ds = psana.DataSource(ds_string)
psana_epix = psana.Detector('ePix100_1')
evr_det = psana.Detector('evr0')

if not args.no_xtcav:
    XTCAVRetrieval=LasingOnCharacterization() 

#grab the first event to get the array sizes
img = psana_epix.photons(next(ds.events()))
n0 = img.shape[0]
n1 = img.shape[1]
nf = args.eventsperfile
imgr = ndimage.rotate(img, args.rotation, order=0)
xes = np.mean(imgr[:,args.rowmin:args.rowmax],axis=1)
ns = xes.shape[0]
total = np.zeros((n0,n1))
total_thresholded = np.zeros((n0,n1))
total_thresholded_l = np.zeros((n1,n0))
nevents_proc = 0
nexth5(args,rank)
# main loop
for ie, evt in enumerate(ds.events()):
    if ie%size!=rank: continue
    if ie%100 == rank:
        print("rank %d processing event %d"%(rank,ie))    
    if ie>=args.max_events:
        break    
    # Time to move on to the next file
    #if dsidx !=0 and (args.eventsperfile is not None and dsidx%args.eventsperfile ==0):
    if dsidx >= args.eventsperfile:
        # Close existing file
        h5sum = h5out.create_dataset('total',(n0,n1), dtype='f4', data=total)
        h5sum_thresh = h5out.create_dataset('total',(n0,n1), dtype='f4', data=total_thresholded)
        #h5sum = total
        
        h5out.close()
        # Increment file name
        subfile += 1
        nexth5(args, rank)
    if evt is None:
        continue
    event_codes = evr_det.eventCodes(evt)
    if event_codes is None:
        continue
    laser_on = 0
    if args.laser_on_code in event_codes:
        laser_on = 1
    #img = psana_epix.image(evt)
    img = psana_epix.photons(evt)
    if img is None:
        continue
    # We have data to store in the frame
   # img = img[args.rowmin:args.rowmax]
   #img = img[args.rowmin:args.rowmax, args.colmin:args.colmax]
    #calculate 1D spectrum
    imgr = ndimage.rotate(img, args.rotation, order=0)
    xes = np.mean(imgr[:,args.rowmin:args.rowmax],axis=1)
    
    # Running sum
    total += img
    tmp = np.copy(img)
    tmp[tmp<args.threshold] = 0
#    total_thresholded += tmp
    #print(rank,idx,dsidx,((size+1)*dsidx)%args.eventsperfile)
    #print(rank, idx)
    #epix[dsidx]=img
    laser[dsidx]=laser_on
    spectra[dsidx]=xes
    total_thresholded_l[dsidx] += tmp
    mtime[dsidx]=evt.get(psana.EventId).time()[0]
    mtimenano[dsidx]=evt.get(psana.EventId).time()[1]
    fiducial[dsidx]=evt.get(psana.EventId).fiducials()
    index[dsidx]=ie
    if not args.no_xtcav:
        if XTCAVRetrieval.processEvent(evt):
            power[dsidx] = np.amax(XTCAVRetrieval.xRayPower()[1])
            agreement[dsidx] = XTCAVRetrieval.reconstructionAgreement()
            try:
                delay= XTCAVRetrieval.pulseFWHM()
            except IndexError:
                delay=None
            if delay is not None:
                pulse[dsidx] = delay[0]
        else:
            power[dsidx]=0
            agreement[dsidx]=0
            pulse[dsidx]=0
    if args.store_raw:
        raw = psana_epix.raw(evt) 
        epixraw[dsidx]=raw
        epix[dsidx]=img
    dsidx += 1
    nevents_proc += 1
print("Rank %d processed %d events"%(rank,nevents_proc))
total_events_proc = comm.reduce(nevents_proc, op=MPI.SUM, root = 0)
if rank == 0:
    print("All ranks processed %d events"%(total_events_proc))
# Close the file if we get this far
total_img = comm.reduce(total, op=MPI.SUM, root=0)
#total_thresholded_img = comm.reduce(total_thresholded, op=MPI.SUM, root=0)
total_thresholded_img = comm.reduce(total_thresholded_l, op=MPI.SUM, root=0)
h5sum = h5out.create_dataset('2D_sum',(n0,n1), dtype='f4',data=total_img)
h5sum_thresh = h5out.create_dataset('2D_sum_thresholded',(n0,n1), dtype='f4',data=total_thresholded_img)
#h5sum = total_img
#h5sum_thresh = total_thresholded_img
h5out.close()
if rank == 0:
    import matplotlib.pyplot as plt
    plt.imshow(total_thresholded_img, vmax=np.median(total_thresholded_img)*2)
    plt.show()
comm.barrier()
if rank == 0:
    vindex = h5py.VirtualLayout(shape=(ie,), dtype='i')
    vlaser = h5py.VirtualLayout(shape=(ie,), dtype='i')
    vspectra = h5py.VirtualLayout(shape=(ie ,ns), dtype='f4')
    vtotal = h5py.VirtualLayout(shape=(ie, ns), dtype='f4')
    vmtime = h5py.VirtualLayout(shape=(ie,), dtype='f4')
    vmtimenano = h5py.VirtualLayout(shape=(ie,), dtype='f4')
    vfiducial = h5py.VirtualLayout(shape=(ie,), dtype='f4')
    if not args.no_xtcav:
        vpower = h5py.VirtualLayout(shape=(ie,), dtype='f4')
        vpulse = h5py.VirtualLayout(shape=(ie,), dtype='f4')
        vagreement = h5py.VirtualLayout(shape=(ie,), dtype='f4')
    if args.store_raw:
        vepix = h5py.VirtualLayout(shape=(ie,n0,n1), dtype='f4')
        vepixraw = h5py.VirtualLayout(shape=(ie,n0,n1), dtype='f4')
    h5filenames=[]
    for h5filename in os.listdir(args.outdir):
        if h5filename.startswith("r%04d-%s" %(args.run, args.tag)) and h5filename.endswith("h5"):
            h5filenames.append(h5filename)
    os.chdir(args.outdir)
    for h5filename in h5filenames:
        print("Processing"+ h5filename)
        h5file = h5py.File(os.path.join(h5filename), 'r')
        indexes = h5file['/index'][:].astype(int)
        for local_index, index in enumerate(indexes):
            if index != 0:
                vindexsource = h5py.VirtualSource(h5file['/index'])
                vindex[index-1] = vindexsource[local_index]
                vlasersource = h5py.VirtualSource(h5file['/laser'])
                vlaser[index-1] = vlasersource[local_index]
                vspectrasource = h5py.VirtualSource(h5file['/spectra'])
                vspectra[index-1] = vspectrasource[local_index]
                vtotalsource = h5py.VirtualSource(h5file['/total_thresholded_l'])
                vtotal[index-1] = vtotalsource[local_index]
                vmtimesource = h5py.VirtualSource(h5file['/machineTime'])
                vmtime[index-1] = vmtimesource[local_index]
                vmtimenanosource = h5py.VirtualSource(h5file['/machineTimeNano'])
                vmtimenano[index-1] = vmtimenanosource[local_index]
                vfiducialsource = h5py.VirtualSource(h5file['/fiducials'])
                vfiducial[index-1] = vfiducialsource[local_index]
                if not args.no_xtcav:
                    vpowersource = h5py.VirtualSource(h5file['/xtcav/power'])
                    vpower[index-1] = vpowersource[local_index]
                    vpulsesource = h5py.VirtualSource(h5file['/xtcav/pulseFWHM'])
                    vpulse[index-1] = vpulsesource[local_index]
                    vagreementsource = h5py.VirtualSource(h5file['/xtcav/agreement'])
                    vagreement[index-1] = vagreementsource[local_index]
                if args.store_raw:
                    vepixsource = h5py.VirtualSource(h5file['/epix'])
                    vepix[index-1] = vepixsource[local_index]
                    vepixrawsource = h5py.VirtualSource(h5file['/epixraw'])
                    vepixraw[index-1] = vepixrawsource[local_index]
        h5file.close()
    vdsfile = h5py.File("VDS-r%04d-%s.h5" %(args.run, args.tag), 'w')
    vdsfile.create_virtual_dataset('index', vindex, fillvalue=0) 
    vdsfile.create_virtual_dataset('laser', vlaser, fillvalue=-1)
    vdsfile.create_virtual_dataset('spectra', vspectra, fillvalue=0)
    vdsfile.create_virtual_dataset('total_thresholded', vtotal, fillvalue=0)
    vdsfile.create_virtual_dataset('machineTime', vmtime, fillvalue=0)
    vdsfile.create_virtual_dataset('machineTimeNano', vmtimenano, fillvalue=0)
    vdsfile.create_virtual_dataset('fiducials', vfiducial, fillvalue=0)
    if not args.no_xtcav:
        vdsfile.create_virtual_dataset('xtcav/power', vpower, fillvalue=0)
        vdsfile.create_virtual_dataset('xtcav/pulseFWHM', vpulse, fillvalue=0)
        vdsfile.create_virtual_dataset('xtcav/agreement', vagreement, fillvalue=0)
    if args.store_raw:
        vdsfile.create_virtual_dataset('epix', vepix, fillvalue=0)
        vdsfile.create_virtual_dataset('epixraw', vepixraw, fillvalue=0)
    vdsfile.close()

