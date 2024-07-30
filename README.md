# Polygon-informed-cross-track-altimetry-PICTA

## About
This repository contains the MATLAB scripts used to generate the PICTA-derived river water level profiles in the work of Ehlers et al. (2024) *Polygon-Informed Cross-Track Altimetry (Picta): Estimating River Water Level Profiles with the Sentinel-6 Altimeter*, see reference below and details in the source publication.

## Getting started
To get a local copy up and running follow these steps.

### Installation
1. Make a local copy of this repository in a directory of your choice, e.g. via 
   ```sh
   git clone https://github.com/...
   ```
2. ```Software/Loadcommonsetting.m``` is used to save and load some info about your local directory. Therein, adjust the ```home_dir``` variable to state the directory of your local copy.
3. Besides built-in MATLAB functions, we also need ```kml2struct``` (https://www.mathworks.com/matlabcentral/fileexchange/35642-kml2struct). Paste the function into the ```Software/``` folder.

### Getting the data
Next we will need to obtain the underlying FFSAR-processed altimetry data from the respective 4TU.ResearchData repository, see https://www.doi.org/10.4121/304db898-f99c-490a-97c4-13f919ae3c05.

3. Unpack the contents of ```Garonne_FFSAR_altimetry_data.zip``` to the directory ```Data/L1b_Garonne/nc/```
4. Unpack the contents of ```Creuse_FFSAR_altimetry_data.zip``` to the directory ```Data/L1b_Creuse/nc/```

## Usage
The PICTA river retracking is implemented here simply as a sequence of 



### References
Ehlers, Frithjof and Slobbe, Cornelis and Schlembach, Florian and Kleinherenbrink, Marcel and Verlaan, Martin, Polygon-Informed Cross-Track Altimetry (Picta): Estimating River Water Level Profiles with the Sentinel-6 Altimeter. Available at SSRN: https://ssrn-com.tudelft.idm.oclc.org/abstract=4851452 or http://dx.doi.org.tudelft.idm.oclc.org/10.2139/ssrn.4851452 

James Slegers (2024). kml2struct (https://www.mathworks.com/matlabcentral/fileexchange/35642-kml2struct), MATLAB Central File Exchange. Retrieved July 30, 2024. 
