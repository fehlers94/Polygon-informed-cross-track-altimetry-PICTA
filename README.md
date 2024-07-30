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
   in the directory ```Data/L2_Garonne/``` to download the appropriate Level-2 altimetry files for the Garonne river, then unzip the files.
7. Execute
   ```sh
   eumdac download -c EO:EUM:DAT:0841 --start 2021-01-01  --end 2022-12-31  --bbox 0.9 46.6 1.0 46.7 --relorbit 73
   ```
   in the directory ```Data/L2_Creuse/``` to and download the appropriate Level-2 altimetry files for the Creuse river, then unzip these files as well.

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
### Step 1.a Generate (uncorrected) river water level profiles from FFSAR data and river polygons 
Open the script ```PICTA_river_retracking.m``` in MATLAB and set the ```river_name``` variable to choose the river scenario to be processed:
```matlab
% choose river scenario to process
river_name = 'Garonne';
%river_name = 'Creuse';
```
Upon execution, the script will iterate over all the FFSAR data in ```Data/L1b_Garonne/nc/``` and save the resulting river water level profiles in ```Results/L2_Garonne/``` as workspace variable (mat-files).

### Step 1.b (optional) Export (uncorrected) river water level profiles to netcdf
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
Upon execution, the script will iterate over all the uncorrected river water level data (mat-files) in ```Results/L2_Garonne/```, apply the geophysical corrections to the river water level estimates and save the data including the corrections to new mat-files in ```Results/L2_Garonne_cor/```.

### Step 3. Export river water level profiles to netcdf
Open the script ```PICTA_export_to_netcdf.m``` in MATLAB and set again the ```river_name``` variable to choose the river scenario to be processed, set ```corrections_included = true```:
```matlab
% choose river scenario
river_name = 'Garonne';
%river_name = 'Creuse';
corrections_included = true;
```
Upon execution, the script will iterate over all the river water level data (mat-files) in ```Results/L2_Garonne_cor/``` and add corresponding netcdf-files with detailed variable descriptions.

## Results
### List of netcdf variables
The netcdf-files produced by ```PICTA_export_to_netcdf.m``` contain the following variables, describing the river water level profiles and some intermediate results of the retracking:
```sh
netcdf file:/home/fehelers/PhD%20Delft/Projects/hydrology%20river/Results/L2_Garonne_cor/S6A_P4_1A_HR______20210108T224244_20210108T233901_20220507T234826_3377_006_070_035_EUM__REP_NT_F06.SEN6_TUDelft.nc
dimensions:
    river_bank_index = 15;
    waveform_index = 55820;
    time = 20;
variables:
 double river_bank_index(river_bank_index=15);
   :description = "Index of river bank. Note that the river bank information is formatted as NaN-separated pairs throughout the 2D arrays, so only field entries with river bank index 1,2 & 4,5 & 7,8 & 10,11 & 13,14 will contain information.";
   :units = "integer";

 double waveform_index(waveform_index=55820);
   :description = "This index represents the satellite positions / waveform number in along track direction. The posting rate of the initial FFSAR waveforms is roughly 1 m.";
   :units = "integer";

 char start_time(time=20);
   :description = "UTC Time of the overpass. Strictly speaking, this is the time corresponding to the waveform index 1, but the satellite records the whole scene within a few seconds anyhow (7 km/s velocity).";
   :units = "yyyy-MM-dd hh:mm:ss";

 double xr(waveform_index=55820, river_bank_index=15);
   :description = "the absolute cross-track distance of the river banks with respect to the local ENU reference frame";
   :units = "meters";

 double along_track(waveform_index=55820, river_bank_index=15);
   :description = "the along-track coordinate of the river banks in the local ENU reference frame (for debugging, should always be ~zero)";
   :units = "meters";

 double cross_track(waveform_index=55820, river_bank_index=15);
   :description = "the cross-track coordinate of the river banks in the local ENU reference frame";
   :units = "meters";

 double East(waveform_index=55820, river_bank_index=15);
   :description = "the eastward coordinate of the river banks with respect to the local ENU reference frame";
   :units = "meters";

 double North(waveform_index=55820, river_bank_index=15);
   :description = "the northward coordinate of the river banks with respect to the local ENU reference frame";
   :units = "meters";

 double Lat(waveform_index=55820, river_bank_index=15);
   :description = "latitude of river banks over WGS84 reference ellipsoid";
   :units = "degrees North";

 double Lon(waveform_index=55820, river_bank_index=15);
   :description = "longitude of river banks over WGS84 reference ellipsoid";
   :units = "degrees East";

 double segmentID(waveform_index=55820, river_bank_index=15);
   :description = "integer number marking which part/segment of the polygon the river banks belong to";
   :units = "integer";

 double yr0(waveform_index=55820, river_bank_index=15);
   :description = "height offset between ground track and river bank location due to Earths curvature in local ENU coordinate system. This is equivalent to xr^2/(2*RE) with Earth radius RE for a round Earth, but evaluated here for the reference WGS84 ellipsoid.";
   :units = "meters";

 double initial_guess_range_gates(waveform_index=55820, river_bank_index=15);
   :description = "initial guess of range gate indices (floating point) based on optimization procedure of initial river values (offsets, slopes)";
   :units = "range gate indices";

 double retracker_range_gates(waveform_index=55820, river_bank_index=15);
   :description = "retracked range gate indices (floating point) using the subwaveforms centered around the initial guess";
   :units = "range gate indices";

 double retracker_river_width(waveform_index=55820, river_bank_index=15);
   :description = "width of retracked river echo (in range gate indices)";
   :units = "range gate indices";

 double expected_river_width(waveform_index=55820, river_bank_index=15);
   :description = "width of initial guess river echo (in range gate indices)";
   :units = "range gate indices";

 double overlap(waveform_index=55820, river_bank_index=15);
   :description = "flag indicating whether one or more subwaveforms at the waveform index are overlapping, which may render the retracked results invalid: 1 = overlap, 0 = no overlap";
   :units = "integer";

 double H_wgs84(waveform_index=55820, river_bank_index=15);
   :description = "measured altitude / elevation of the individual river banks over WGS84 (averaging over opposite banks NOT applied, geophysical corrections NOT applied)";
   :units = "meters";

 double yr_retracked_ENU(waveform_index=55820, river_bank_index=15);
   :description = "intermediate height measurement (meters) in local ENU coordinate system: yr_retracked_ENU + yr0 = H_wgs84";
   :units = "meters";

 double R_retracked(waveform_index=55820, river_bank_index=15);
   :description = "retracked slant range between satellite and river bank (averaging over opposite banks NOT applied, geophysical corrections NOT applied)";
   :units = "meters";

 double H_wgs84_avg(waveform_index=55820, river_bank_index=15);
   :description = "measured altitude / elevation of the individual river banks over WGS84 (averaging over opposite banks applied, no geophysical corrections applied)";
   :units = "meters";

 double Lat_avg(waveform_index=55820, river_bank_index=15);
   :description = "The latitude of the river segments. The river segments are actually lines along the cross-track direction, of which this coordinate represents only the center.";
   :units = "degrees North";

 double Lon_avg(waveform_index=55820, river_bank_index=15);
   :description = "The longitude of the river segments. The river segments are actually lines along the cross-track direction, of which this coordinate represents only the center.";
   :units = "degrees East";

 double subwf_len(waveform_index=55820, river_bank_index=15);
   :description = "length of subwaveform in units of range gate indices";
   :units = "range gate indices";

 double mask_valid(waveform_index=55820, river_bank_index=15);
   :description = "A preliminary quality flag based on the retracking results. A value of 1 marks potentially valid measurements, a value of 0 marks potentially invalid measurements.";
   :units = "integer";

 double H_wgs84_cor(waveform_index=55820, river_bank_index=15);
   :description = "measured altitude / elevation of the individual river banks over WGS84 (averaging over opposite banks NOT applied, geophysical corrections applied)";
   :units = "meters";

 double H_wgs84_avg_cor(waveform_index=55820, river_bank_index=15);
   :description = "measured altitude / elevation of the individual river banks over WGS84 (averaging over opposite banks applied, geophysical corrections applied)";
   :units = "meters";

 double Lat_sat(waveform_index=55820);
   :description = "The latitude of the satellite CoM (WGS84)";
   :units = "degrees North";

 double Lon_sat(waveform_index=55820);
   :description = "The longitude of the satellite CoM (WGS84)";
   :units = "degrees East";

 double Alt_sat(waveform_index=55820);
   :description = "The altitude of the satellite CoM above the WGS84 reference ellipsoid.";
   :units = "meters";

 double tracker_range(waveform_index=55820);
   :description = "The calibrated tracker range measured from the CoM of the satellite platform, as in the EUMETSAT L1a files, but interpolated to FFSAR posting rate.";
   :units = "meters";
```
### Visualization of the data
![image](https://github.com/user-attachments/assets/74bebef8-71ca-475f-ae30-8d9cfcbf39b8)

The matlab script producing this visualization of a single river water level profile is provided below.
```matlab
filename = '/home/fehelers/PhD Delft/Projects/hydrology river/Results/L2_Garonne_cor/S6A_P4_1A_HR______20210527T182205_20210527T191822_20220509T142357_3377_020_070_035_EUM__REP_NT_F06.SEN6_TUDelft.nc'
water_level = ncread(filename,"H_wgs84_avg_cor");
lon = ncread(filename,"Lon_avg");
lat = ncread(filename,"Lat_avg");
lon_sat = ncread(filename,"Lon_sat");
lat_sat = ncread(filename,"Lat_sat");
mask = ncread(filename,"mask_valid");
mask = (mask==1);
% ##### plot the map with elevations of the river segments #####
fig = figure('units','inch','position',[0,0,8,8],'visible','on');
set(gcf,'color','w');

geoplot(lat_sat,lon_sat,'LineWidth',2,'LineStyle','--','Color','k'); hold on
geoscatter(lat(:),lon(:),'k.');hold on
geoscatter(lat(mask),lon(mask),[],water_level(mask),'filled');hold on

geolimits([44.3 44.6],[0.05 0.16])
caxis([50,70])
geotickformat('dd')
colorbar()
colormap(rgb)
legend('satellite track','river points','river height')
```

### References
Ehlers, Frithjof and Slobbe, Cornelis and Schlembach, Florian and Kleinherenbrink, Marcel and Verlaan, Martin, Polygon-Informed Cross-Track Altimetry (Picta): Estimating River Water Level Profiles with the Sentinel-6 Altimeter. Available at SSRN: https://ssrn-com.tudelft.idm.oclc.org/abstract=4851452 or http://dx.doi.org.tudelft.idm.oclc.org/10.2139/ssrn.4851452 

James Slegers (2024). kml2struct (https://www.mathworks.com/matlabcentral/fileexchange/35642-kml2struct), MATLAB Central File Exchange. Retrieved July 30, 2024. 
