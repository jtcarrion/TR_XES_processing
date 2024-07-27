import h5py
import sys
import numpy as np
from psana import *


parser = argparse.ArgumentParser(description='Script for idenitifying and filtering events based on laser activation. Output laser on/off fiducials in H5 datasets')
parser.add_argument('-e', dest='exp', type=str, required=True, help='experiment string (i.e. mfxlw0953)')
parser.add_argument('-r', dest='run', type=int, required=True, help='run number to process')
args = parser.parse_args()

experiment_id = args.exp
run_number = args.run
ds = DataSource('exp=%s:run=%d'%(experiment_id, run_number)+':smd')
evr_det = Detector('evr0')
laser_on_code = 210
laser_on_fid = np.zeros(1)
laser_off_code = 211
laser_off_fid = np.zeros(1)
on_cnt = 0
off_cnt = 0

for event_number, event in enumerate(ds.events()):
    ec = evr_det.eventCodes(event)
    if ec is None:
        continue
    laser_status = None

    eid = event.get(EventId)
    fid = eid.fiducials()
    sec = eid.time()[0]
    nsec = eid.time()[1]

    if laser_on_code in ec:
        print(event_number, 'Laser on\tF: ', fid, sec, nsec)
        laser_status = True
        laser_on_fid = np.append(laser_on_fid, fid)
#        on_cnt += 1
    if laser_off_code in ec:
        print(event_number, 'Laser off\tF: ', fid, sec, nsec)
        laser_status = False
        laser_off_fid = np.append(laser_off_fid, fid)
#        off_cnt += 1
    if laser_status is None:
        print(event_number, 'Laser status unknown')

on_fid=laser_on_fid[laser_on_fid != 0]
off_fid=laser_off_fid[laser_off_fid != 0]
print(len(on_fid), len(off_fid))

h5 = h5py.File('laser_on_fiducials.h5', 'w')
h5.create_dataset('laser_on_fiducials', data=on_fid) 
h5.close()

h5_off = h5py.File('laser_off_fiducials.h5', 'w')
h5_off.create_dataset('laser_off_fiducials', data=off_fid)
h5_off.close()
