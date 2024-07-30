# Polygon-informed-cross-track-altimetry-PICTA

## About
This repository contains the MATLAB scripts used to generate the PICTA-derived river water level profiles in the work of Ehlers et al. (2024) *Polygon-Informed Cross-Track Altimetry (Picta): Estimating River Water Level Profiles with the Sentinel-6 Altimeter*, see reference below and details in the source publication.

## Getting started
In order to run these scripts you will need an installation of MATLAB (R2021a) including the mapping toolbox.
To get a local copy of this repository up and running to reproduce the results of the paper, follow these steps.

### Installation
1. Make a local copy of this repository in a directory of your choice, e.g. via 
   ```sh
   git clone https://github.com/...
   ```
2. The script ```Software/Loadcommonsetting.m``` is used to save and load some info about your local directory. Therein, adjust the ```home_dir``` variable to state the directory of your local repository, e.g.
   ```matlab
   home_dir = '/home/fehelers/PhD Delft/Projects/Polygon-informed-cross-track-altimetry-PICTA';
   ```
4. Besides built-in MATLAB functions, we also need ```kml2struct``` (https://www.mathworks.com/matlabcentral/fileexchange/35642-kml2struct). Download the function into the ```Software``` folder.

### Getting all the data
Next we will need to obtain the underlying FFSAR-processed altimetry data over the Garonne and Creuse rivers from the respective 4TU.ResearchData repository, see https://www.doi.org/10.4121/304db898-f99c-490a-97c4-13f919ae3c05.

4. Unpack the contents of ```Garonne_FFSAR_altimetry_data.zip``` to the directory ```Data/L1b_Garonne/nc/```
5. Unpack the contents of ```Creuse_FFSAR_altimetry_data.zip``` to the directory ```Data/L1b_Creuse/nc/```

Additionally, we will have to obtain the corresponding Level-2 altimetry data from EUMETSAT (http://doi.org/10.15770/EUM_SEC_CLM_0096), which contain the values of the geophysical corrections and which are not contained in the lower Level-1a and Level-1b products. The data can be obtained via the EUMETSAT Data Access Client (EUMDAC) client, see https://user.eumetsat.int/resources/user-guides/eumetsat-data-access-client-eumdac-guide on how to set it up. Once you have set up EUMDAC (made a user account, set your credentials, etc.) you can 

6. execute
   ```sh
   eumdac download -c EO:EUM:DAT:0841 --start 2021-01-01  --end 2022-12-31  --bbox -1 42 1 50 --relorbit 35  
   ```
   in the directory ```Data/L2_Garonne/``` to download the appropriate Level-2 altimetry files for the Garonne river.
7. Execute
   ```sh
   eumdac download -c EO:EUM:DAT:0841 --start 2021-01-01  --end 2022-12-31  --bbox 0.9 46.6 1.0 46.7 --relorbit 73
   ```
   in the directory ```Data/L2_Creuse/``` to repeat the same step for the Creuse river.


## Brief software description
The PICTA river retracking is implemented here simply as a sequence of three rudimentary MATLAB scripts: ```PICTA_river_retracking.m```, ```PICTA_apply_geophysical_corrections.m``` and```PICTA_export_to_netcdf.m```, the headers of which contain detailed information about their purpose, inputs, outputs and variable descriptions.
```matlab
% =========================================================================
% Script Name: PICTA_river_retracking.m
% -------------------------------------------------------------------------
% Purpose: 
%   This script reads FFSAR-processed altimetry data and a polygon of the
%   Garonne or Creuse rivers in France and applies the PICTA method to 
%   derive dense river water level profiles (~10 m along-track resolution).
% ...
```

```matlab
=========================================================================
% Script Name: PICTA_apply_geophysical_corrections.m
% -------------------------------------------------------------------------
% Purpose: 
%   This script reads PICTA-processed river water level profile data produced by
%   'PICTA_river_retracking.m' and applies the geophysical corrections interpolated from
%   the corresponding EUMETSAT Level-2 files, see
%       "EUMETSAT for Copernicus (2023): Poseidon-4 Altimetry Level 2 High 
%       Resolution (baseline version F08) - Sentinel-6 - Reprocessed, 
%       European Organisation for the Exploitation of Meteorological Satellites, 
%       DOI: 10.15770/EUM_SEC_CLM_0096,
%       https://dx.doi.org/10.15770/EUM_SEC_CLM_0096"```
% ...
```

```matlab
% =========================================================================
% Script Name: PICTA_export_to_netcdf.m
% -------------------------------------------------------------------------
% Purpose: 
%   This script reads PICTA-processed river water level profile data produced by
%   "PICTA_river_retracking.m" or
%   "PICTA_apply_geophysical_corrections.m"
%   and exports the data to netcdf.
% ...
```


## Usage
### Step 1.a Generate river water level profiles from FFSAR data and river polygons 
Open the script ```PICTA_river_retracking.m``` in MATLAB and set the ```river_name``` variable to choose the river scenario to be processed:
```matlab
% choose river scenario to process
river_name = 'Garonne';
%river_name = 'Creuse';
```
Upon execution, the script will iterate over all the FFSAR data in ```Data/L1b_Garonne/nc/``` and save the resulting river water level profiles in ```Results/L2_Garonne/``` as workspace variable (mat-files).

### Step 1.b (optional) Export river water level profiles to netcdf
Open the script ```PICTA_export_to_netcdf.m``` in MATLAB and set again the ```river_name``` variable to choose the river scenario to be processed, keep ```corrections_included = false```:
```matlab
% choose river scenario
river_name = 'Garonne';
%river_name = 'Creuse';
corrections_included = false;
```
Upon execution, the script will iterate over all the river water level data (mat-files) in ```Results/L2_Garonne/``` and add corresponding netcdf-files with detailed variable descriptions.

### Step 2. Apply geophysical corrections 
Open the script ```PICTA_apply_geophysical_corrections.m``` in MATLAB and set again the ```river_name``` variable to choose the river scenario to be processed:
```matlab
% choose river scenario
river_name = 'Garonne';
%river_name = 'Creuse';
```
Upon execution, the script will iterate over all the raw river water level data (mat-files) in ```Results/L2_Garonne/```, apply the geophysical corrections to the river water level estimates and save the data including the corrections to new mat-files in ```Results/L2_Garonne_cor/```.

### Step 3. Export river water level profiles to netcdf
Open the script ```PICTA_export_to_netcdf.m``` in MATLAB and set again the ```river_name``` variable to choose the river scenario to be processed, set ```corrections_included = true```:
```matlab
% choose river scenario
river_name = 'Garonne';
%river_name = 'Creuse';
corrections_included = true;
```
Upon execution, the script will iterate over all the river water level data (mat-files) in ```Results/L2_Garonne_cor/``` and add corresponding netcdf-files with detailed variable descriptions.

### References
Ehlers, Frithjof and Slobbe, Cornelis and Schlembach, Florian and Kleinherenbrink, Marcel and Verlaan, Martin, Polygon-Informed Cross-Track Altimetry (Picta): Estimating River Water Level Profiles with the Sentinel-6 Altimeter. Available at SSRN: https://ssrn-com.tudelft.idm.oclc.org/abstract=4851452 or http://dx.doi.org.tudelft.idm.oclc.org/10.2139/ssrn.4851452 

James Slegers (2024). kml2struct (https://www.mathworks.com/matlabcentral/fileexchange/35642-kml2struct), MATLAB Central File Exchange. Retrieved July 30, 2024. 
