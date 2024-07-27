# XES_data_vizualization
This script can process, sort, and visualize X-ray Emission Spectroscopy (XES) data collected from SLAC National Accelerator Laboratory. Pairing XES and SFX data allows for observations of quantum changes in protein complexes at a femtosecond time scale. This script is meant for processing on SLAC's HPC utilizing Psana and Python 3.9. Raw data collected from the Epix100 XES and XTCAV SFX dector in the LCLS MFX hutch are processed and stored as Virtural Data Set (VDS) structures and the post-processing is stored in H5 datasets. 

## Installation:
```sh
git clone https://github.com/jtcarrion/TR_XES_processing.git
```

## Usage:
`mpiXESVDSandXTCAV.py`: Main script for processing raw Epix100 XES data from LCLS MFX   
`filter_laser_events.py`: Filtering events based on laser on/off  
`filter_cheetah_hits.py`: Filtering events based on SFX crystal hits  
`plot_xes.py`: Visualize and quantify the XES spectra 

The workflow above can processe and sort serial-femtosecond XES data based on laser On/Off events as well as crystal hits/no hits, with femtosecond percision. With the current state of the detector, analyzing laser on vs off events was only possible before sorting hits vs non hits to reduce the background noice enough and see the correct peaks. Pairing with the `filter_cheetah_hits.py` script, a significant increase is seen in the signal to noise ratio, allowing for more percise quantification. 
