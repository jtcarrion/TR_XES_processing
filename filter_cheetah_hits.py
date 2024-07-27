import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse


parser = argparse.ArgumentParser(description='Script for processing H5 datasets with combined Cheetah Xtal Hits')
parser.add_argument("-f","--files", nargs="+", help="xes h5 file list, space separated")
parser.add_argument("-w","--window", default=1, type=int, help="window size for smoothing (running mean, default=1, no smoothing)")
parser.add_argument("-o","--output", type=str, default=None, help="output filename for png plot")
parser.add_argument("-c","--CheetahHits", nargs="+", help="lst file containing hits from Cheetah")
args = parser.parse_args()

# Create a numpy array of the Xtal hit fiducials from Cheeta, stored as int
hits = open("%s" %args.cheetaHits[0]).read()
print('----------')
hits_f = hits.split()
hits_fiducials = np.array(hits_f, dtype=str)

#print('-------------------')
#print(len(hits_fiducials))

nf = len(args.files)
laser = []
spectra = []

# Loop through various VDS files to average over several runs
for i in range(nf):
    h5 = h5py.File(args.files[i])
    spectra.append(h5['spectra'][:])
    laser.append(h5['spectra'][:])

laser = h5['laser'][:]

fiducials = h5['fiducials'][:]
machineTimeNano = h5['machineTimeNano'][:]
machineTime = h5['machineTime'][:]

print('-------------------')
print("shape of timestamp from h5")
print(fiducials.shape)
fiducials_array = np.array(fiducials, dtype= int)
print('-------------------')

# Remove the unnecesary information from the entire timestamp
hits_fid = [x.split('-')[2] for x in hits_fiducials]
hits_mtn = [x.split('-')[1] for x in hits_fiducials]
hits_mt =  [x.split('-')[0] for x in hits_fiducials]

print('-------------------')
print("shape of fiducials after split")
print(len(hits_fid))
hits_fid = [int(x) for x in hits_fid] # Store fiducials as int
print(type(hits_fid[0]))
print('-------------------')

spectra = h5['spectra'][:]
spectra = np.flip(spectra)
print(spectra.shape[0], 'total frames')

# Get the index where the laser is On vs Off
won = np.where(laser == 0)
woff = np.where(laser == 1)

print('-------------------')
print("won length")
print(len(won[0]))
print('-------------------')
print("woff length")
print(len(woff[0]))

hits_index = []
no_hits_index = []

# Loop through each fiducial from the VDS and if it is a hit, append the INDEX of that event to their respective array
z = 0
for x in fiducials_array:
    if x in hits_fid:
        hits_index.append(z)
    else:
        no_hits_index.append(z)
    z+=1

print('-------------------')
print("Number of hit indexes")
print(len(hits_index))
print('-------------------')

won_hits = []
woff_hits = []

# Loop through each hit event/index and sort based on laser On vs Off
for y in hits_index:
    if y in won[0]: won_hits.append(y)
    if y in woff[0]: woff_hits.append(y)

# In Experimnet LY0022 and LW9515, should see 67% Laser On vs 33% Laser Off
print('-----------')
print('----len(won_hits)-------')
print(len(won_hits))
print('----len(woff_hits)-------')
print(len(woff_hits))
print('-----------')

# Create new spectra where the index is equivalent to only Laser On Hits, Laser Off Hits, All Hits, and no hits
spectra_on_hits = np.squeeze(spectra[won_hits, :])
spectra_off_hits = np.squeeze(spectra[woff_hits, :])
spectra_hits = np.squeeze(spectra[hits_index, :])
spectra_no_hits = np.squeeze(spectra[no_hits_index, :])
n1 = spectra_on_hits.shape[0]
n2 = spectra_off_hits.shape[0]
n3 = spectra_hits.shape[0]
n4 = spectra_no_hits.shape[0]

#print('-------------------')
#print("n1 is :" + str(n1))
#print("n2 is :" + str(n2))
#print("n3 is :" + str(n3))
#print("n4 is :" + str(n4))

# Get the average spectra over these events mentioned above
# Sum up the spectra and then divide by the number of events in that respective group
spectra_on_hits_mean = np.sum(spectra_on_hits, axis=0)/n1
np.savetxt('XES_spectra_on_hits_mean.txt', spectra_on_hits_mean)

spectra_off_hits_mean = np.sum(spectra_off_hits, axis=0)/n2
np.savetxt('XES_spectra_off_hits_mean.txt', spectra_off_hits_mean)

spectra_hits_mean = np.sum(spectra_hits, axis=0)/n3
np.savetxt('XES_spectra_hits_mean.txt', spectra_hits)

spectra_no_hits_mean = np.sum(spectra_no_hits, axis=0)/n4
np.savetxt('XES_spectra_no_hits_mean.txt', spectra_no_hits)

# Use non Xtal hits as background correction
spectra_hits_diff = np.subtract(spectra_hits_mean, spectra_no_hits_mean)

print(spectra_on_hits.shape[0], 'laser on frames')
print(spectra_off_hits.shape[0], 'laser off frames')
print(spectra_hits.shape[0], 'Hit frames')
print(spectra_no_hits.shape[0], 'No hit frames')

# Plot the resulting Spectra for each group
# Fig 1. Laser On vs Off
plt.figure(1)
plt.plot(spectra_on_hits_mean, label='laser on')
plt.plot(spectra_off_hits_mean, label='laser off')
plt.xlabel('Intensity')
plt.ylabel('pixel')
plt.title('On vs Off for crytals hits')
#plt.plot(np.mean(spectra_off_mean), label='laser off')
plt.legend()
plt.savefig('XES_spectra_on_off.png', dpi=150)
plt.show()

# Fig 2. No Xtal hit spectra
plt.figure(2)
plt.plot(spectra_no_hits_mean, label='no hits')
#plt.imshow(spectra)
plt.xlabel('Intensity')
plt.ylabel('pixel')
plt.title('No hits mean')
#np.savetxt('xes-1Dspec.txt', xes_aft_avg_norm)
plt.savefig('XES_spectra_no_hits.png', dpi=150)
plt.show()

# Fig 3. Xtal hits spectra
plt.figure(3)
#plt.plot(laser, '.')
plt.plot(spectra_hits_mean, label='hits')
#plt.show()
plt.xlabel('Intensity')
plt.ylabel('pixel')
plt.title('Hits mean')
#np.savetxt('xes-1Dspec.txt', xes_aft_avg_norm)
#plt.savefig('XES_spectra_on_and_off_hits.png' , dpi=150)
plt.show()

# Fig 4. Xtal hits minus the no hit/background
plt.figure(4)
plt.plot(spectra_hits_diff, label='Hits vs background')
#plt.show()
plt.xlabel('Intensity')
plt.ylabel('pixel')
plt.title('Hits spectra minus no hits spectra')
#np.savetxt('xes-1Dspec.txt', xes_aft_avg_norm)
plt.savefig('XES_spectra_on_and_off_hits.png' , dpi=150)
plt.show()

