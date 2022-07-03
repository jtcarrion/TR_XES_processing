# XES_data_vizualization
This script can process, sort, and visualize X-ray Emission Spectroscopy (XES) data collected from SLAC National Accelerator Laboratory. Pairing XES and serial 
femtosecond protein crystallography allows us to observe quantum changes in protein complexes at a femtosecond time scale. 

This script is meant for processing on SLAC's HPC utilizing Psana and Python 3. This script has only been tested with data collected from the Epix100 XES dector in the  LCLS MFX hutch. 

This script takes in one or many Virtual Datasets (VDS) from the previous processing script as well as a list of crystal hits from Cheeta. The VDS will contain the XES data collected from the dector and the list of hits will be provided by a 3rd party hit finder program like Cheeta

The script here will sort the processed XES data based on laser On vs Off events as well as crystal hits vs no hit events with femtosecond percision. The Epix detector had odd background noise during data collection and the injection media for the crystals were new so vizualizing the extact events during a crystal hit will allow for limited background noise as well as allow us to subtract the non crystal hits data for further background subtraction. 

With the current state of the detector, analyzing laser on vs off events was only possible after sorting hits vs non hits to reduce the background noice enough and see the correct peaks. 
