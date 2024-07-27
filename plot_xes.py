import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Script for visualizing XES spectras')
parser.add_argument('-f', dest='data', type=str, required=True, help='Input H5 file')
parser.add_argument('--upper_bound', dest='upper', type=int, required=True, help='Upper bound for the XES spectra')
parser.add_argument('--lower_bound', dest='lower', type=int, required=True, help='Lower bound for the XES spectra')
args = parser.parse_args()

img = None
files = args.data
nf = len(files)
for i in range(nf):
    print(files[i])
    h5 = h5py.File(files[i])
    if img is None:
        img = h5['2D_sum_thresholded'][:]
    else:
        img += h5['2D_sum_thresholded'][:]


xes = np.sum(img[:,args.lower_bound:args.upper_bound],axis=1)

vmin = np.percentile(img.ravel(),10)
vmax = np.percentile(img.ravel(),99.5)
plt.imshow(img, vmin=vmin,vmax=vmax)
plt.savefig('2D_sum.png',dpi=150)
plt.show()

plt.plot(xes)
plt.savefig('1D_xes.png',dpi=150)
plt.show()
np.savetxt('1D_xes.txt', xes)
