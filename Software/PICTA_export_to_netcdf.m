% =========================================================================
% Script Name: PICTA_export_to_netcdf.m
% -------------------------------------------------------------------------
% Purpose: 
%   This script reads PICTA-processed river water level profile data produced by
%   "PICTA_apply_geophysical_corrections.m" and exports the data to netcdf
%
% Usage:
%   Run the script directly in the MATLAB environment.
%
% Inputs:
%   The script assumes that the intermediate .mat-files containing the PICTA-
%   processed Sentinel-6 data including geophysical corrections (output of "PICTA_apply_geophysical_corrections.m") 
%   is located in the folders
%   - Results/L2_Creuse_cor
%   or
%   - Results/L2_Garonne_cor
%
% Outputs:
%   The script produces netcdf files with identical names as the input files 
%   and adds them to the directories
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
%   This code is licensed under the MIT License.
% =========================================================================

clear all
LoadCommonSettings 

% choose river scenario
%river_name = 'Creuse';
river_name = 'Garonne';
corrections_included = false;

if corrections_included
    fnames = dir(fullfile(PathRSLT,['L2_' river_name '_cor'],'*.mat'))
else
    fnames = dir(fullfile(PathRSLT,['L2_' river_name],'*.mat'))
end

var_names_2d_units_description = [
    {'xr'                       } {'meters'             } {'the absolute cross-track distance of the river banks with respect to the local ENU reference frame'};
    {'along_track'              } {'meters'             } {'the along-track coordinate of the river banks in the local ENU reference frame (for debugging, should always be ~zero)'};
    {'cross_track'              } {'meters'             } {'the cross-track coordinate of the river banks in the local ENU reference frame'};
    {'East'                     } {'meters'             } {'the eastward coordinate of the river banks with respect to the local ENU reference frame'};
    {'North'                    } {'meters'             } {'the northward coordinate of the river banks with respect to the local ENU reference frame'};
    {'Lat'                      } {'degrees North'      } {'latitude of river banks over WGS84 reference ellipsoid'};
    {'Lon'                      } {'degrees East'       } {'longitude of river banks over WGS84 reference ellipsoid'};
    {'segmentID'                } {'integer'            } {'integer number marking which part/segment of the polygon the river banks belong to'};
    {'yr0'                      } {'meters'             } {'height offset between ground track and river bank location due to Earths curvature in local ENU coordinate system. This is equivalent to xr^2/(2*RE) with Earth radius RE for a round Earth, but evaluated here for the reference WGS84 ellipsoid.'};
    {'initial_guess_range_gates'} {'range gate indices' } {'initial guess of range gate indices (floating point) based on optimization procedure of initial river values (offsets, slopes)'};
    {'retracker_range_gates'    } {'range gate indices' } {'retracked range gate indices (floating point) using the subwaveforms centered around the initial guess'};
    {'retracker_river_width'    } {'range gate indices' } {'width of retracked river echo (in range gate indices)'};
    {'expected_river_width'     } {'range gate indices' } {'width of initial guess river echo (in range gate indices)'};
    {'overlap'                  } {'integer'            } {'flag indicating whether one or more subwaveforms at the waveform index are overlapping, which may render the retracked results invalid: 1 = overlap, 0 = no overlap'};
    {'H_wgs84'                  } {'meters'             } {'measured altitude / elevation of the individual river banks over WGS84 (averaging over opposite banks NOT applied, geophysical corrections NOT applied)'};
    {'yr_retracked_ENU'         } {'meters'             } {'intermediate height measurement (meters) in local ENU coordinate system: yr_retracked_ENU + yr0 = H_wgs84'};
    {'R_retracked'              } {'meters'             } {'retracked slant range between satellite and river bank (averaging over opposite banks NOT applied, geophysical corrections NOT applied)'};
    {'H_wgs84_avg'              } {'meters'             } {'measured altitude / elevation of the individual river banks over WGS84 (averaging over opposite banks applied, no geophysical corrections applied)'};
    {'Lat_avg'                  } {'degrees North'      } {'The latitude of the river segments. The river segments are actually lines along the cross-track direction, of which this coordinate represents only the center.'};
    {'Lon_avg'                  } {'degrees East'       } {'The longitude of the river segments. The river segments are actually lines along the cross-track direction, of which this coordinate represents only the center.'};
    {'subwf_len'                } {'range gate indices' } {'length of subwaveform in units of range gate indices'};
    {'mask_valid'               } {'integer'            } {'A preliminary quality flag based on the retracking results. A value of 1 marks potentially valid measurements, a value of 0 marks potentially invalid measurements.'};
    {'H_wgs84_cor'              } {'meters'             } {'measured altitude / elevation of the individual river banks over WGS84 (averaging over opposite banks NOT applied, geophysical corrections applied)'};
    {'H_wgs84_avg_cor'          } {'meters'             } {'measured altitude / elevation of the individual river banks over WGS84 (averaging over opposite banks applied, geophysical corrections applied)'}]

%{'waveform_index'           } {'integer'            } {'This index represents the satellite positions / waveform number in along track direction. The posting rate of the initial FFSAR waveforms is roughly 1 m.'};

var_names_1d_units_description = [
    {'Lat_sat'                  } {'degrees North'      } {'The latitude of the satellite CoM (WGS84)'};
    {'Lon_sat'                  } {'degrees East'       } {'The longitude of the satellite CoM (WGS84)'};
    {'Alt_sat'                  } {'meters'             } {'The altitude of the satellite CoM above the WGS84 reference ellipsoid.'};
    {'tracker_range'            } {'meters'             } {'The calibrated tracker range measured from the CoM of the satellite platform, as in the EUMETSAT L1a files, but interpolated to FFSAR posting rate.'}]

for i=1:numel(fnames)
    fname = fnames(i).name
    load(fullfile(fnames(i).folder,fname))
    

    res.mask_valid = 1*res.mask_valid;
    res.waveform_index = res.wf_idx;
    res.Lat_sat = res.GEO.LAT';
    res.Lon_sat = res.GEO.LON';
    res.Alt_sat = res.GEO.H';
    res.tracker_range = res.MEA.tracker_range';

    % Define the dimensions
    nx = size(res.Lat, 1);
    ny = size(res.Lat, 2);

    % Create a NetCDF file
    %ncid = netcdf.create(fullfile(PathRSLT,['L2_' river_name '_cor'],[fname(1:end-4) '_TUDelft.nc']), 'NC_WRITE');
    ncid = netcdf.create(fullfile(fnames(i).folder,[fname(1:end-4) '.nc']), 'NC_WRITE');

    %%
    % Define dimensions
    dimid_x = netcdf.defDim(ncid, 'river bank index', nx);
    dimid_y = netcdf.defDim(ncid, 'waveform index', ny);
    dimid_time = netcdf.defDim(ncid, 'time', 20); % Define time dimension

    varid_x = netcdf.defVar(ncid, 'river bank index', 'double', dimid_x);
    varid_y = netcdf.defVar(ncid, 'waveform index', 'double', dimid_y);
    varid_time = netcdf.defVar(ncid, 'start time', 'char', dimid_time);

    netcdf.putAtt(ncid, varid_time, 'description', 'UTC Time of the overpass. Strictly speaking, this is the time corresponding to the waveform index 1, but the satellite records the whole scene within a few seconds anyhow (7 km/s velocity).');
    netcdf.putAtt(ncid, varid_time, 'units', 'yyyy-MM-dd hh:mm:ss');

    netcdf.putAtt(ncid, varid_x, 'description', 'Index of river bank. Note that the river bank information is formatted as NaN-separated pairs throughout the 2D arrays, so only field entries with river bank index 1,2 & 4,5 & 7,8 & 10,11 & 13,14 will contain information.');
    netcdf.putAtt(ncid, varid_x, 'units', 'integer');

    netcdf.putAtt(ncid, varid_y, 'description', 'This index represents the satellite positions / waveform number in along track direction. The posting rate of the initial FFSAR waveforms is roughly 1 m.');
    netcdf.putAtt(ncid, varid_y, 'units', 'integer');

    % Define all other variables with descriptions
    varids_2d = [];
    for k = 1:size(var_names_2d_units_description,1)
        varids_2d.(var_names_2d_units_description{k,1}) = netcdf.defVar(ncid, var_names_2d_units_description{k,1}, 'double', [dimid_x, dimid_y]);
        netcdf.putAtt(ncid, varids_2d.(var_names_2d_units_description{k,1}), 'description', var_names_2d_units_description{k,3});
        netcdf.putAtt(ncid, varids_2d.(var_names_2d_units_description{k,1}), 'units', var_names_2d_units_description{k,2});
    end

    varids_1d = [];
    for k = 1:size(var_names_1d_units_description,1)
        varids_1d.(var_names_1d_units_description{k,1}) = netcdf.defVar(ncid, var_names_1d_units_description{k,1}, 'double', [dimid_y]);
        netcdf.putAtt(ncid, varids_1d.(var_names_1d_units_description{k,1}), 'description', var_names_1d_units_description{k,3});
        netcdf.putAtt(ncid, varids_1d.(var_names_1d_units_description{k,1}), 'units', var_names_1d_units_description{k,2});
    end

    netcdf.endDef(ncid);

    %% 
    % Write data to variables

    netcdf.putVar(ncid, varid_x, 1:nx);
    netcdf.putVar(ncid, varid_y, 1:ny);
    netcdf.putVar(ncid, varid_time, char(res.GEO.Start_Time)); % Write time value

    for k = 1:size(var_names_2d_units_description,1)
        if corrections_included
            netcdf.putVar(ncid, varids_2d.(var_names_2d_units_description{k,1}), res.(var_names_2d_units_description{k,1}));
        else
            if ~contains(var_names_2d_units_description{k,1},'_cor')
                netcdf.putVar(ncid, varids_2d.(var_names_2d_units_description{k,1}), res.(var_names_2d_units_description{k,1}));
            end
        end
    end

    for k = 1:size(var_names_1d_units_description,1)
        netcdf.putVar(ncid, varids_1d.(var_names_1d_units_description{k,1}), res.(var_names_1d_units_description{k,1}));
    end

    % Close the NetCDF file
    netcdf.close(ncid);
end

% %% some test-reading and plotting
% %clear all
% 
% % wse = ncread(['nc/' fname(1:end-4) '_TUDelft.nc'],"water_surface_elevation")
% % lon = ncread(['nc/' fname(1:end-4) '_TUDelft.nc'],"longitude")
% % lat = ncread(['nc/' fname(1:end-4) '_TUDelft.nc'],"latitude")
% 
% 
% wse = ncread(fullfile(fnames(1).folder,[fname(1:end-4) '.nc']),"H_wgs84_avg");
% lon = ncread(fullfile(fnames(i).folder,[fname(1:end-4) '.nc']),"Lon_avg");
% lat = ncread(fullfile(fnames(i).folder,[fname(1:end-4) '.nc']),"Lat_avg");
% lon_sat = ncread(fullfile(fnames(i).folder,[fname(1:end-4) '.nc']),"Lon_sat");
% lat_sat = ncread(fullfile(fnames(i).folder,[fname(1:end-4) '.nc']),"Lat_sat");
% mask = ncread(fullfile(fnames(i).folder,[fname(1:end-4) '.nc']),"mask_valid");
% 
% scatter(lon(mask==1),lat(mask==1),[],wse(mask==1),'filled');hold on
% plot(lon_sat,lat_sat)
