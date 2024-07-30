% =========================================================================
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
%       https://dx.doi.org/10.15770/EUM_SEC_CLM_0096"
%
% Usage:
%   Run the script directly in the MATLAB environment.
%
% Inputs:
%   The script will compare the file-identifiers of PICTA and EUMETSAT
%   Level 2 data for satellite cycle and relative orbit number. The script assumes that 
%   the corresponding EUMETSAT Level 2-files are downloaded and unpacked into the
%   directories
%   - Data/L2_Creuse
%   - Data/L2_Garonne
%   and that the PICTA-processed data is readily located in
%   - Results/L2_Creuse
%   - Results/L2_Garonne
%   after execution of "PICTA_river_retracking.m"
%
% Outputs:
%   The output is similar to that of "PICTA_river_retracking.m", namely a struct array 
%   'res' describing the river water level profile. However, this script
%   adds two new fields, "H_wgs84_cor" and "H_wgs84_avg_cor", where the
%   subscript "_cor" indicates that the geophysical corrections have been
%   applied.
%   By default, the data is saved as a intermediate .mat file into the folders
%   - Results/L2_Creuse_cor
%   or
%   - Results/L2_Garonne_cor
%
% Author: 
%   Frithjof Ehlers
%   Department of Geoscience and Remote Sensing
%   Faculty of Civil Engineering
%   Delft University of Technology
%   f.ehlers@tudelft.nl
%
% Date of Creation:
%   July 11, 2024
%
% License:
%   This code is licensed under the XXX License.
% =========================================================================

clear all
LoadCommonSettings

%river_name = 'Creuse'
river_name = 'Garonne'

DOM = [-60 60; -180 180]

% directory of own processing and eumetsat data containing correction values
L2_dir_own = fullfile(PathRSLT,['L2_' river_name]);
L2_dir_eum = fullfile(PathDATA,['L2_' river_name]);

L2_files_own   = dir([L2_dir_own '/S6*.mat']);
%%
for i = 1:numel(L2_files_own)
    L2_name_own = L2_files_own(i).name;
    id = L2_name_own(71:end-18);
    
    % find matching L2 file
    L2_match   = dir([L2_dir_eum '/*' id '*/*' id '*/S6A_P4_2__HR_STD*.nc']);
    L2_match = fullfile(L2_match.folder,L2_match.name);
    
    % load EUMETSAT L2 file corrections for inland water
    [~,COR] = S6_read_L2(L2_match,DOM);
    
    % load own L2 data struct
    load(fullfile(L2_dir_own,L2_name_own))
    
    % correct entries in res
    
    if isfield(res,'GEO')
        x = res.GEO.LAT';
    else
        x = res.Lat_avg;
    end
    
    total_cor = COR.model_dry_tropo_cor_measurement_altitude(x) + ...   
                COR.model_wet_tropo_cor_measurement_altitude(x) + ...   
                COR.iono_cor_gim(x) + ...
                COR.solid_earth_tide(x) + ...
                COR.pole_tide(x) + ...% geophysical corrections according to DINARDO DISS
                COR.model_instr_cor_range_ocean(x); % since this modelled range cor was missing in L1a tracker range
    
    %res.H_wgs84_avg([3,6,9,12,15],:) = NaN;
    %res.H_wgs84([3,6,9,12,15],:) = NaN;
            
    res.H_wgs84_cor = res.H_wgs84 - total_cor % added in range but substracted in altitude
    res.H_wgs84_avg_cor = res.H_wgs84_avg - total_cor % added in range but substracted in altitude
    
    % save the res as the same L2 struct
    save(fullfile([L2_dir_own '_cor'],L2_name_own),'res','-v7.3')
end
