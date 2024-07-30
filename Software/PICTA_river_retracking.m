% =========================================================================
% Script Name: PICTA_river_retracking.m
% -------------------------------------------------------------------------
% Purpose: 
%   This script reads FFSAR-processed altimetry data and a polygon of the
%   Garonne or Creuse rivers in France and applies the PICTA method to 
%   derive dense river water level profiles (~10 m along-track resolution).
%
% Usage:
%   Run the script directly in the MATLAB environment.
%
% Inputs:
%   The script will use the following .kml-files containing the river 
%   polygon information:
%   - Data/OSM/garonne.kml
%   - Data/OSM/la_creuse.kml
%   Additionally, the FFSAR-processed Sentinel-6 altimeter data in form of
%   .nc-files should be located in 
%   - Data/L1b_Garonne/nc
%   - Data/L1b_Creuse/nc
%   An example of a file name of the FFSAR-processed data is
%   "S6A_P4_1A_HR______20210121T195334_20210121T204951_20230518T112742_3377_007_146_073_EUM__REP_NT_F08.SEN6_TUDelft.nc",
%   which is essentially a copy of the file-identifiers from the underlying Level 1a
%   EUMETSAT data collection "EO:EUM:DAT:0838", see
%       "EUMETSAT for Copernicus (2023): Poseidon-4 Altimetry Level 1A High
%       Resolution (baseline version F08) - Sentinel-6 - Reprocessed, European 
%       Organisation for the Exploitation of Meteorological Satellites, DOI: 
%       10.15770/EUM_SEC_CLM_0093, https://dx.doi.org/10.15770/EUM_SEC_CLM_0093"
%
% Outputs:
%   The output is a struct array 'res' for each of the overpasses with the 
%   data fields describing the river water level profile listed below. Be
%   aware that no geophysical corrections are applied to the range/elevation at this step.
%   By default, the data is saved as a intermediate .mat file into the folders
%   - Results/L2_Creuse
%   or
%   - Results/L2_Garonne
%   The first dimension of each field denotes the river bank index, the second dimension
%   marks the waveform index. For better data handling, the river banks are saved
%   in NaN-separated pairs, since two banks belong to a single segment. I.e. 
%   if there are two river segments (s1 and s2) inside one waveform i, then, e.g., 
%   the longitude info of both their banks (b1 and b2) is ordered as
%   Lon(:,i) = [s1_b1 s1_b2 NaN s2_b1 s2_b2 NaN NaN NaN NaN ...]
%   This ordering is consistent throughout the 2D array.
%   We implicitly assume here that no more than 5 river segments are in
%   view at any time.
%                            xr: [15×55820 double] - the absolute cross-track distance of the river banks (meters)
%                   along_track: [15×55820 double] - the along-track distance of the river banks with respect to the waveform index (meters, for debugging, should always be ~zero)
%                   cross_track: [15×55820 double] - the cross-track distance of the river banks (meters)
%                          East: [15×55820 double] - river bank position equivalent to along-track and cross-track distance, but rotated to local East-North(-Up) coordinate system (meters)
%                         North: [15×55820 double] - river bank position equivalent to along-track and cross-track distance, but rotated to local East-North(-Up) coordinate system (meters)
%                           Lat: [15×55820 double] - Latitude of river banks over WGS84 reference ellipsoid (degrees North)
%                           Lon: [15×55820 double] - Longitude of river banks over WGS84 reference ellipsoid (degrees East)
%                     segmentID: [15×55820 double] - integer number marking which larger polygon segment the river banks belong to
%                           yr0: [15×55820 double] - height offset between
%                           ground track and river bank location due to
%                           Earths curvature in local ENU coordinate system
%                           (meters). This is equivalent to xr^2/(2*RE)
%                           with Earth radius RE for a round Earth, but
%                           evaluated for the reference WGS84 ellipsoid
%     initial_guess_range_gates: [15×55820 double] - initial guess of range gate indices (floating point) based on optimization procedure of initial river values (offsets, slopes)
%         retracker_range_gates: [15×55820 double] - retracked range gate indices (floating point) using the subwaveforms centered around the initial guess
%         retracker_river_width: [15×55820 double] - width of retracked river echo (in range gate indices)
%          expected_river_width: [15×55820 double] - width of initial guess river echo (in range gate indices)
%                       overlap: [15×55820 double] - flag indicating whether one or more subwaveforms at the waveform index are overlapping
%                       H_wgs84: [15×55820 double] - measured altitude / elevation of the individual river banks over WGS84 (meters; no averaging applied, no geophysical corrections applied)
%              yr_retracked_ENU: [15×55820 double] - intermediate height measurement (meters) in local ENU coordinate system: yr_retracked_ENU + yr0 = H_wgs84
%                   R_retracked: [15×55820 double] - retracked slant range between satellite and river bank (meters; no averaging applied, no geophysical corrections applied)
%                   H_wgs84_avg: [15×55820 double] - measured altitude / elevation of the river segments over WGS84 (meters; averaged over opposite river banks, no geophysical corrections applied)
%                       Lat_avg: [15×55820 double] - Latitudes of river segment centers over WGS84 reference ellipsoid (degrees North)
%                       Lon_avg: [15×55820 double] - Longitudes of river segment centers over WGS84 reference ellipsoid (degrees East)
%                     subwf_sum: [15×55820 double] - summed intensity over subwaveform
%                     subwf_len: [15×55820 double] - length of subwaveform in units of range gate indices
%                        wf_idx: [15×55820 double] - along-track waveform index (integer)
%                           GEO: [1×1 struct]      - struct containing the satellites altitude/longitude/latitude data of the satellite
%                           MEA: [1×1 struct]      - struct containing the satellites tracker_range setting
%                    mask_valid: [15×55820 logical]- preliminary validity flag: 1=valid, 0=invalid
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
warning('off','all')

% choose river scenario to process
river_name = 'Garonne';
%river_name = 'Creuse';

% initialize the range resolution and refbin
% some S6 parameters
dr = 0.18974206202531646; % range gate spacing with 512 range gates (2x zero-padding)
refBin = 257; % reference range gate within waveform (1-based indexing)

di = 1; % spacing for along-track operations in units of waveform index

% set some specific river scenario information
if strcmp(river_name,'Garonne')
    contains_string = '_070_035_EUM'; % case-specific file ending of EUMETSAT L1A files
    
    % Read the KML file with river polygons into a struct:
    river_segments = kml2struct(fullfile(PathDATA,'OSM','garonne.kml'));
    % some preprocessing to throw away the smaller rivers:
    river_segments([3 4 5 6 7 8 9 10])=[];
    
    % define segments, each segment will get its own offset/slope for fitting
    segments=[];
    segments(1).LonLat = [-0.5 44; -0.5,44.8 ; 0.5,44.8; 0.5,44.; -0.5 44];
    
    % range gate buffer for choosing the sub-waveforms
    buffer = 20;
    
    seg_thr = 15; % any cross-track river segment narrower than 15 m is omitted
    
elseif strcmp(river_name,'Creuse')
    contains_string = '_146_073_EUM'; % case-specific file ending of EUMETSAT L1A files
    
    % Read the KML file with river polygons into a struct:
    river_segments = kml2struct(fullfile(PathDATA,'OSM','la_creuse.kml'));
    
    % define segments, each segment will get its own offset/slope for fitting
    segments=[];
    segments(1).LonLat = [0.83 46.785; 0.91,46.742; 0.915,46.718; 1.016,46.617; 1.31 46.56; 1.31 46.8; 0.83 46.8; 0.83 46.785];
    segments(2).LonLat = [0.83 46.785; 0.91,46.742; 0.915,46.718; 1.016,46.617; 1.31 46.56; 0.4 46.56; 0.83 46.785];
    segments(3).LonLat = [0.83 46.785; 0.83 47.2; 0.1 47.2;  0.1 46.785; 0.83 46.785;];

    % range gate buffer for choosing the sub-waveforms
    buffer = 20;
    
    seg_thr = 15; % any cross-track river segment narrower than 15 m is omitted
end

%% collect all polygon segments in polyshape vector and merge them (take union)
for i = 1:numel(river_segments)
    if i ==1
        polyvec = [polyshape(river_segments(i).Lat,river_segments(i).Lon)];
    end
    % polyvec(i) = 
    polyvec(i) = polyshape(river_segments(i).Lat,river_segments(i).Lon);
end

% merge polygon into one
polyout = union(polyvec);
% remove the holes within the river polygon (e.g. small islands within the river)
polyout = rmholes(polyout);

% collect latitude+longitude in a struct called 'river'
river.Lat = polyout.Vertices(:,1);
river.Lon = polyout.Vertices(:,2);

% get FFSAR file information from data directory
files   = dir([PathDATA,'/L1b_' river_name '/nc/*_P4_1A_HR' '*' contains_string '*.nc']);

%% PICTA processing iterated over all files in directory
tic
for file_num = 1:numel(files)
    toc
    filename = files(file_num).name
    
    %% load S6 FFSAR data
    fprintf('### loading FFSAR data ### \n')
    filepath = fullfile(files(file_num).folder,filename);
    CS1b = [];
    % for variable descriptions see netcdf file
    CS1b.SAR.data = ncread(filepath,"waveform_power");
    CS1b.GEO.LON = ncread(filepath,"longitude");
    CS1b.GEO.LAT = ncread(filepath,"latitude");
    CS1b.GEO.H = ncread(filepath,"altitude");
    CS1b.MEA.tracker_range = ncread(filepath,"tracker_range");
    CS1b.GEO.Start_Time = datetime(ncread(filepath,"start_time")');
    
    CS1b.SAR.data10 = conv2(CS1b.SAR.data,ones(1,10),'same'); % moving average of 10 waveforms ~ 10 meters along-track
    CS1b.SAR.data20 = conv2(CS1b.SAR.data,ones(1,20),'same'); % moving average of 20 waveforms ~ 20 meters along-track

    %% find the intersection points of cross-track FFSAR footprints and river polygon (Lat,Lon)
    fprintf('### finding river intersection points ###  \n')
    
    % initialize struct array of intersection points with the river
    river_points = [];
    Nmax = numel(CS1b.GEO.LAT);
    river_points(Nmax).xr = []; % positve cross-track distance of intersection point (positive)
    river_points(Nmax).along_track = []; % along-track coordinate relative to waveform ground point (must be ~zero, for testing)
    river_points(Nmax).cross_track = []; % cross-track coordinate relative to waveform ground point (with sign)
    river_points(Nmax).East = []; % equivalent to along_track, cross_track, but in North-East coordinate system
    river_points(Nmax).North = []; % equivalent to along_track, cross_track, but in North-East coordinate system
    river_points(Nmax).Lat = []; % latitudes of intersection points
    river_points(Nmax).Lon = []; % longitudes of intersection points
    river_points(Nmax).segmentID = []; % segmentID=1,2,3... defines the segment to which the intersection point belongs
    river_points(Nmax).yr0 = []; % height offset of intersection points due to ellipsoids curvature (replacing to x^2/(2*RE), see paper)
    
    %tic
    
    % identify whether river is a multipolygon polygon by identifying
    % the different regions described by a single polygon
    multiple_regions_flag = (polyshape(river.Lat,river.Lon).NumRegions>1);

    if multiple_regions_flag
        idxs = find(isnan(river.Lat));
        for i = 1:numel(idxs)+1
            if i==1
                region(i).start = 1;
                region(i).end = idxs(i)-1;
            elseif i==numel(idxs)+1
                region(i).start = region(i-1).end+2; % skipping the NaN value itself
                region(i).end = numel(river.Lat);
            else 
                region(i).start = region(i-1).end+2; % skipping the NaN value itself
                region(i).end = idxs(i)-1;
            end
        end
    end

    % iterate over all waveforms (ground track locations) to find
    % cross-track intersection points with river polygon
    for i = 1:di:numel(CS1b.GEO.LAT)-1
        %toc
        %p = i/numel(CS1b.GEO.LAT)*100

        % Transform into ENU coordinates
        wgs84 = wgs84Ellipsoid;
        lat0 = CS1b.GEO.LAT(i);
        lon0 = CS1b.GEO.LON(i);
        h0 = 0;

        [riverE,riverN,riverZ] = geodetic2enu(river.Lat,river.Lon,0,lat0,lon0,h0,wgs84);
        [atrackE,atrackN,trackZ] = geodetic2enu(CS1b.GEO.LAT(i+1),CS1b.GEO.LON(i+1),0,lat0,lon0,h0,wgs84);
        ctrackE = atrackN;
        ctrackN = -atrackE;

        % normalize atrack and ctrack vectors
        scale = 1;
        ctrackNn = scale*ctrackN./sqrt(ctrackN^2+ctrackE^2);
        ctrackEn = scale*ctrackE./sqrt(ctrackN^2+ctrackE^2);
        atrackNn = scale*atrackN./sqrt(atrackN^2+atrackE^2);
        atrackEn = scale*atrackE./sqrt(atrackN^2+atrackE^2);

        %turn the coordinate system into A and C coordinates (along track, cross track)
        if strcmp(river_name,'Garonne')
            theta = acos(atrackEn);
        elseif strcmp(river_name,'Lek')
            theta = -acos(atrackEn);
        elseif strcmp(river_name,'Creuse')
            theta = acos(atrackEn);
        end

        [riverA,riverC] = rotate2d(riverE,riverN,theta);

        if ~multiple_regions_flag
            % mask the polygon to the relevant snippet (assuming that the polygon is sampled denser than 500 m)
            % this deletion of points from the polygon causes troubles in case of multiple regions
            mask = (riverA < 1000)&((riverA > -1000));

            % find intercection points (slow)
            [in,out] = intersect(polyshape(riverA(mask),riverC(mask)),[-0 -12e3; 0 12e3]); % considering a -12,+12 km FFSAR footprint (value not important as long as it is chosen large enough to cover the range window)
        else
            in = [];
            %figure;
            for k = 1:numel(region)
                vecA = riverA(region(k).start:region(k).end);
                vecC = riverC(region(k).start:region(k).end);
                % mask the polygon to the relevant snippet (assuming that the polygon is sampled denser than 500 m)
                mask = (vecA < 1000)&((vecA > -1000));
                
                if sum(mask)>4 % require more than four points to consider the region at all
                    % find intercection points (slow)
                    [in_region,out_region] = intersect(polyshape(vecA(mask),vecC(mask)),[-0 -12e3; 0 12e3]); % considering a -12,+12 km FFSAR footprint (value not important as long as it is chosen large enough to cover the range window)
                    if ~isempty(in_region)
                        if isempty(in)
                            in = in_region;
                        else
                            in = [in; NaN NaN; in_region];
                        end
                    end
                end
            end
        end

        if ~isempty(in)
            if isempty(river_points(i).xr)
                % fill the arrays with the obtained values
                river_points(i).along_track = in(:,1);
                river_points(i).cross_track = in(:,2);
                river_points(i).xr = sqrt(river_points(i).cross_track.^2 + river_points(i).along_track.^2);

                % redo the turn of variables
                [river_points(i).East,river_points(i).North] = rotate2d(river_points(i).along_track,river_points(i).cross_track,-theta);
                [river_points(i).Lat,river_points(i).Lon,river_points(i).yr0] = enu2geodetic(river_points(i).East,river_points(i).North,0,lat0,lon0,h0,wgs84);

                river_points(i).segmentID = zeros(size(river_points(i).xr));
                for seg_idx = 1:numel(segments)
                    river_points(i).segmentID = river_points(i).segmentID + seg_idx*inpolygon(river_points(i).Lat,river_points(i).Lon,segments(seg_idx).LonLat(:,2),segments(seg_idx).LonLat(:,1));
                end
            else
                % fill the arrays with the obtained values, but separate
                % them with a NaN from the earlier obtained values
                along_track = cat(1,river_points(i).along_track,[NaN;in(:,1)]);
                cross_track = cat(1,river_points(i).cross_track,[NaN;in(:,2)]);
                river_points(i).xr = cat(1,river_points(i).xr,[NaN; sqrt(cross_track.^2 + along_track.^2)]);

                % redo the turn of variables
                [East,North] = rotate2d(along_track,cross_track,-theta);
                %river_points(i).East = cat(1,river_points(i).East,[NaN;East]);
                %river_points(i).North = cat(1,river_points(i).North,[NaN;North]);

                [Lat,Lon,yr0] = enu2geodetic(East,North,0,lat0,lon0,h0,wgs84);
                river_points(i).Lat = cat(1,river_points(i).Lat,[NaN;Lat]);
                river_points(i).Lon = cat(1,river_points(i).Lon,[NaN;Lon]);
                river_points(i).yr0 = cat(1,river_points(i).yr0,[NaN;yr0]);
            end
        end
        
    end

    %% objective function optimization to find optimal height and slope of each river segment over wgs84 prior to retracking
    fprintf('### optimizing first guess of river elevation ### \n')
    
    if strcmp(river_name,'Garonne')
        % define initial guess of offset and slopes over lat and lon for each segment
        H0 = 55;
        dH_lat=-35;
        dH_lon=25;
        x0 = [H0 dH_lat dH_lon];
        
        % define objective function (see end of this script) and parameters
        fun = @(x)ObjFunc2(x,river_points,CS1b,river,numel(segments),[]);
        options = optimset('Display','iter','PlotFcns',@optimplotfval);

        %[x,fval,exitflag,output] = fminsearch(fun,x0,options)
        [x,fval,exitflag,output] = fminsearch(fun,x0);
        
    elseif strcmp(river_name,'Creuse')
        x0 = [102 -55 45 102 -55 45 100 -50 42];
        equality_constraint=[1 4]; % expressing that x(1) = x(4), the offset of segment 1 and 2 should be equal, as the two segments occupy the same waveforms and since the rivers meet
        
        % define objective function (see end of this script) and parameters
        fun = @(x)ObjFunc2(x,river_points,CS1b,river,numel(segments),equality_constraint);
        options = optimset('Display','iter','PlotFcns',@optimplotfval);

        %[x,fval,exitflag,output] = fminsearch(fun,x0,options)
        [x,fval,exitflag,output] = fminsearch(fun,x0);
        x(equality_constraint(1)) = x(equality_constraint(2));
    end

    
    %% plot the forward model with the optimized parameters:
    
    % convolute the data just for visualization and the cost function
    % efficiency:
    SARdata = CS1b.SAR.data10;
    %SARdata = CS1b.SAR.data_pseudoDD;
%     figure;
%     %x=x0
%     %ax1 = subplot(1,2,1)
%     %imagesc(log10(CS1b.SAR.data_pseudoDD));hold on
%     %imagesc(log10(SARdata));hold on
%     %ax2 = subplot(1,2,2)
%     %imagesc(log10(CS1b.SAR.data_pseudoDD));hold on
%     imagesc(log10(SARdata));hold on
%     %H0 = 43; 
%     %dH = 0;
% %     dH = -0.27/1000;
% %     H0 = 0.5*(H0_Tonneins + H0_Reole)%68; 
%     % dH = 0;
%     %d2H = 0;
% %     H0 = 59.;
% %     dH_lat=-35;
% %     dH_lon=25;
% %     
%     for i = 1:di:numel(river_points)%numel(river_points)
%         % assume we have a river height that varies linearly over wgs84:
%         %H_wgs84 = H0 + dH*river_points(i).Dist;%; + d2H*river_points(i).Dist.^2;
%         %H_wgs84 = H0 + dH_lat*(river_points(i).Lat-mean(river.Lat,'omitnan'))  + dH_lon*(river_points(i).Lon-mean(river.Lon,'omitnan'));
%         H_wgs84 = zeros(size(river_points(i).Lat));
%         for seg_idx = 1:numel(segments)
%             H_wgs84 = H_wgs84 + (seg_idx==river_points(i).segmentID).*(x(3*(seg_idx-1)+1) + x(3*(seg_idx-1)+2)*(river_points(i).Lat-mean(river.Lat,'omitnan'))  + x(3*(seg_idx-1)+3)*(river_points(i).Lon-mean(river.Lon,'omitnan')));
%         end
%         % now we have to substract the elevation that is due to the cross-track
%         % distance over a curved surface, and can then use this height to
%         % calculate the range in the radargram
%         yr = H_wgs84 - river_points(i).yr0;
%         Hs = CS1b.GEO.H(i);
%         xr = river_points(i).xr;
%     %     R0 = ALT + xr.^2/(2*ALT) - CS1b.MEA.tracker_range(i);
%     %     RH = H.*( 1 + d.^2/(2*ALT^2));
%     %     R = R0+RH;
%         R = Hs + xr.^2/(2*Hs) - CS1b.MEA.tracker_range(i) - yr;
%         rg = R/dr+257;
%         plot(i*ones(size(R)),rg,'k.','MarkerSize',3.5);hold on
%         %plot(i*ones(size(R)),rg,'k-');hold on
%         %plot(i*ones(size(R)),rg,'ro');hold on
%         %geoplot(river_points(i).Lat,river_points(i).Lon,'b.');hold on
%     end
%     %linkaxes([ax1 ax2],'xy')

    
    
    %% remove all segments that are smaller than seg_thr
    for i = 1:numel(river_points)
        river_points(i).xr;
        num_seg = (numel(river_points(i).xr)+1)/3;
        for k=1:num_seg
            wc = abs(river_points(i).xr(3*k-2)-river_points(i).xr(3*k-1));
            if wc<seg_thr
                river_points(i).xr(3*k-2) = NaN;
                river_points(i).xr(3*k-1) = NaN;
            end
        end
    end
    
    %% apply threshold retracking and obtain water level over WGS84 from the FFSAR data, using the information in river_points
    fprintf('### performing threshold retracking ### \n')
    for i = 1:di:numel(river_points)
        if ~isempty(river_points(i).xr)
            % compute the expected range gate limits of the river echo (from optimization results x):
            Hinit_wgs84 = zeros(size(river_points(i).Lat));
            for seg_idx = 1:numel(segments)
                Hinit_wgs84 = Hinit_wgs84 + (seg_idx==river_points(i).segmentID).*(x(3*(seg_idx-1)+1) + x(3*(seg_idx-1)+2)*(river_points(i).Lat-mean(river.Lat,'omitnan'))  + x(3*(seg_idx-1)+3)*(river_points(i).Lon-mean(river.Lon,'omitnan')));
            end
            
            yr_init = Hinit_wgs84 - river_points(i).yr0;
            Hs = CS1b.GEO.H(i);
            xr = river_points(i).xr;
            R = Hs + xr.^2/(2*Hs) - CS1b.MEA.tracker_range(i) - yr_init;
            rg = R/dr+257;
            rg_decimal = R/dr+257;

            river_points(i).initial_guess_range_gates = rg;

            rg = round(rg);

%             if strcmp(river_name,'Garonne')
%                 buffer = 10;
%             elseif strcmp(river_name,'Creuse')
%                 buffer = 20;
%             end
            
            % now retrack the sar echos over subwaveforms centered around the expected echo:
            threshold = 0.1;
            wf_domains = [];
            river_points(i).retracker_range_gates = [];

            %for loop over cross-sections
            for k = 1:(numel(rg)+1)/3
                if rg(2*k-2+k)>rg(2*k-1+k)
                    rgs =  round(rg(2*k-1+k)-buffer:rg(2*k-2+k)+buffer);
                else
                    rgs =  round(rg(2*k-2+k)-buffer:rg(2*k-1+k)+buffer);
                end
                rgs(rgs>512)=[];
                rgs(rgs<1)=[];
                
                if ~isnan(rgs)
                    subwf = SARdata(rgs,i);
                    subwf = subwf./max(subwf);
                    wf_domains(k).rgs = rgs;

                    mask = subwf > threshold;
                    idxs = find(mask);

                    % find the peak with highest intensity in case there are multiple peaks
                    split_idxs = find(diff(idxs)>1);
                    if ~isempty(split_idxs)
                        % split into three arrays
                        sub_idxs = [];
                        num_peaks = numel(split_idxs)+1;
                        intensity = zeros(1,num_peaks);

                        for peak = 1:num_peaks
                            if peak == 1
                                sub_idxs(peak).idxs = idxs(1:split_idxs(peak));
                            elseif peak == num_peaks
                                sub_idxs(peak).idxs = idxs(split_idxs(peak-1)+1:end);
                            else
                                sub_idxs(peak).idxs = idxs(split_idxs(peak-1)+1:split_idxs(peak));
                            end
                            % remove peak if it is at the subwaveform border (contains indices 1 or numel(subwf))
                            if ismember(1,sub_idxs(peak).idxs)|ismember(numel(subwf),sub_idxs(peak).idxs)
                                sub_idxs(peak).idxs = [];
                                intensity(peak) = NaN;
                            else
                                intensity(peak) = sum(subwf(sub_idxs(peak).idxs))./numel(sub_idxs(peak).idxs);
                            end
                        end
                        % replace idxs with sub_idxs of highest intensity
                        [~,max_peak] = max(intensity);
                        idxs = sub_idxs(max_peak).idxs;    
                    end

                    % perform threshold retracking with those new idxs
                    if ~isempty(idxs)    
                        ind_l = idxs(1);
                        ind_r = idxs(end);

                        if (ind_l ~= 1) & (ind_r ~= numel(subwf))
                            rg_l = (threshold - subwf(ind_l-1))/(subwf(ind_l)-subwf(ind_l-1)) + rgs(ind_l-1);
                            rg_r = (threshold - subwf(ind_r))/(subwf(ind_r+1)-subwf(ind_r)) + rgs(ind_r);
                        else
                            rg_l = NaN;
                            rg_r = NaN;
                        end
                    else
                        rg_l = NaN;
                        rg_r = NaN;
                    end

                    % save the retracked range gates in array similar to the rest in river_points
                    if rg(2*k-2+k)>rg(2*k-1+k)
                        river_points(i).retracker_range_gates(2*k-2+k) = rg_r;
                        river_points(i).retracker_range_gates(2*k-1+k) = rg_l;
                    else
                        river_points(i).retracker_range_gates(2*k-2+k) = rg_l;
                        river_points(i).retracker_range_gates(2*k-1+k) = rg_r;
                    end

                    % save the retracked river width (in range gate units):
                    river_points(i).retracker_river_width(2*k-2+k) = abs(rg_r-rg_l);
                    river_points(i).retracker_river_width(2*k-1+k) = abs(rg_r-rg_l);
                    % and the expected according to forward model:
                    river_points(i).expected_river_width(2*k-2+k) = abs(rg_decimal(2*k-2+k)-rg_decimal(2*k-1+k));
                    river_points(i).expected_river_width(2*k-1+k) = abs(rg_decimal(2*k-2+k)-rg_decimal(2*k-1+k));
                    
                    % save the sum and length of the subwf
                    river_points(i).subwf_sum(2*k-2+k) = sum(subwf);
                    river_points(i).subwf_sum(2*k-1+k) = sum(subwf);
                    river_points(i).subwf_len(2*k-2+k) = numel(subwf);
                    river_points(i).subwf_len(2*k-1+k) = numel(subwf);
                else
                    river_points(i).retracker_range_gates(2*k-2+k) = NaN;
                    river_points(i).retracker_range_gates(2*k-1+k) = NaN;
                    
                    river_points(i).retracker_river_width(2*k-2+k) = NaN;
                    river_points(i).retracker_river_width(2*k-1+k) = NaN;
                    
                    % and the expected according to forward model:
                    river_points(i).expected_river_width(2*k-2+k) = abs(rg_decimal(2*k-2+k)-rg_decimal(2*k-1+k));
                    river_points(i).expected_river_width(2*k-1+k) = abs(rg_decimal(2*k-2+k)-rg_decimal(2*k-1+k));
                end
                % adjust NaNs and structure
                river_points(i).retracker_range_gates = river_points(i).retracker_range_gates(:);
                river_points(i).retracker_range_gates(river_points(i).retracker_range_gates==0) = NaN;  

                river_points(i).retracker_river_width = river_points(i).retracker_river_width(:);
                river_points(i).retracker_river_width(river_points(i).retracker_river_width==0) = NaN;  

                river_points(i).expected_river_width = river_points(i).expected_river_width(:);
                river_points(i).expected_river_width(river_points(i).expected_river_width==0) = NaN;  
            end

            %detect overlap between waveform domains and set the flag true
            n = numel(wf_domains);
            river_points(i).overlap = false;
            for p = 1:n-1
                for q = p+1:n
                    if ~isempty(intersect(wf_domains(p).rgs,wf_domains(q).rgs))
                        river_points(i).overlap = true;
                    end
                end
            end

            % now reverse those following steps in order to go from range gate to Hwgs84:
            %     yr = H_wgs84 - river_points(i).yr0;
            %     Hs = CS1b.GEO.H(i);
            %     xr = river_points(i).xr;
            %     R = Hs + xr.^2/(2*Hs) - CS1b.MEA.tracker_range(i) - yr;
            %     rg = R/dr+257;

            R_re = dr*(river_points(i).retracker_range_gates - refBin);R_re=R_re(:);
            yr_re = -(R_re - Hs - xr.^2/(2*Hs) + CS1b.MEA.tracker_range(i));yr_re=yr_re(:);
            river_points(i).H_wgs84 = yr_re + river_points(i).yr0;
            river_points(i).yr_retracked_ENU = yr_re;
            river_points(i).R_retracked = R_re;

            % now after the raw calculation, average the heights of a single
            % cross-section and along-river distance:
            for k = 1:(numel(rg)+1)/3
                river_points(i).H_wgs84_avg(2*k-2+k) = 0.5*( river_points(i).H_wgs84(2*k-2+k)+river_points(i).H_wgs84(2*k-1+k) );
                river_points(i).H_wgs84_avg(2*k-1+k) = 0.5*( river_points(i).H_wgs84(2*k-2+k)+river_points(i).H_wgs84(2*k-1+k) );
                
                river_points(i).Lat_avg(2*k-2+k) = 0.5*( river_points(i).Lat(2*k-2+k)+river_points(i).Lat(2*k-1+k) );
                river_points(i).Lat_avg(2*k-1+k) = 0.5*( river_points(i).Lat(2*k-2+k)+river_points(i).Lat(2*k-1+k) );

                river_points(i).Lon_avg(2*k-2+k) = 0.5*( river_points(i).Lon(2*k-2+k)+river_points(i).Lon(2*k-1+k) );
                river_points(i).Lon_avg(2*k-1+k) = 0.5*( river_points(i).Lon(2*k-2+k)+river_points(i).Lon(2*k-1+k) );
            end
            river_points(i).H_wgs84_avg = river_points(i).H_wgs84_avg(:);
            
            river_points(i).Lon_avg = river_points(i).Lon_avg(:);
            river_points(i).Lat_avg = river_points(i).Lat_avg(:);
        end
    end


    %% after retracking, bring the data into a unified format (2D array) for faster plotting and analysis
    % use all variables besides "overlap"
    fprintf('### saving data to /Results/ ### \n')
    
    var_names = fieldnames(river_points);

    max_ct_i = 15; % assuming not more than five river segments in one cross-track footprint

    % initialize new struct
    res = [];
    for v = 1:numel(var_names)
        var = var_names{v};
        res.(var) = NaN*zeros(max_ct_i,numel(river_points));
    end

    % fill the new struct with the respective entries
    for i = 1:numel(river_points)
        for v = 1:numel(var_names)
            var = var_names{v};
            if ~strcmp(var,'overlap')
                max_ind = numel(river_points(i).(var));
                res.(var)(1:max_ind,i) = river_points(i).(var);
            else
                max_ind = numel(river_points(i).('Lat'));
                if max_ind == 0
                    res.(var)(1:max_ind,i) = river_points(i).(var);
                else
                    res.(var)(1:max_ct_i,i) = river_points(i).(var);
                end
            end
        end
    end

    res.wf_idx = (1:numel(river_points)).*ones(max_ct_i,numel(river_points));
    res.GEO = CS1b.GEO;
    res.MEA = CS1b.MEA;
    mask_valid = ~((res.overlap==1)|(abs(res.retracker_river_width - res.expected_river_width - 3) > 4));%|isnan(res.retracker_river_width)
    res.mask_valid = mask_valid;
    
    save(fullfile(PathRSLT,['L2_' river_name],[filename(1:end-3) '.mat']),'res','-v7.3' )
    
    %% plot the retracker range gates above the waveform power:
    fig = figure('units','inch','position',[0,0,18,12],'visible','off');
    set(gcf,'color','w');

    ax1 = subplot(2,1,1)
    imagesc(log10(SARdata));hold on

    ax2 = subplot(2,1,2)
    imagesc(log10(SARdata));hold on

    plot(res.wf_idx(:),res.initial_guess_range_gates(:),'k-','LineWidth',1);hold on
    plot(res.wf_idx(:),res.retracker_range_gates(:),'r-','LineWidth',1);hold on
    plot(res.wf_idx(mask_valid),res.retracker_range_gates(mask_valid),'b-','LineWidth',1);hold on

    legend('valid','invalid')
    colormap(bone)
    
    linkaxes([ax1,ax2],'xy')

    saveas(gcf,fullfile(PathRSLT, ['wf_data_fit_' filename '.png']))
    close all

    %% plot map of averaged river elevations
    fig = figure('units','inch','position',[0,0,12,12],'visible','off');
    set(gcf,'color','w');

    rgb = [ ...
        94    79   162
        50   136   189
       102   194   165
       171   221   164
       230   245   152
       255   255   191
       254   224   139
       253   174    97
       244   109    67
       213    62    79
       158     1    66  ] / 255;

    rgb = interp1(1:11,rgb,1:0.05:11);
    % ############### plot the map with elevations of the river segments

    geoplot(CS1b.GEO.LAT,CS1b.GEO.LON,'LineWidth',2,'LineStyle','--','Color','k'); hold on
    geoplot(river.Lat,river.Lon,'LineWidth',1,'Color','#0072BD'); hold on

    geoscatter(res.Lat_avg(mask_valid),res.Lon_avg(mask_valid),[],res.H_wgs84_avg(mask_valid),'filled');hold on
    %geoscatter(res.Lat_avg(:),res.Lon_avg(:),[],res.H_wgs84_avg(:),'filled');hold on

    if strcmp(river_name,'Garonne')
        geolimits([44.3 44.6],[0.05 0.16])
        caxis([50,70])
    elseif strcmp(river_name,'Creuse')
        geolimits([46.55 47.08],[0.5 1.2])
        caxis([75.,130.])
    end
    geotickformat('dd')
    %geobasemap colorterrain
    geobasemap none
    %colormap(gca,winter)
    colorbar()
    colormap(rgb)
    legend('satellite track','river polygon','river height')
    set(gcf,'color','w');

    saveas(gcf,fullfile(PathRSLT, ['elevations_' filename '.png']))
    close all
end

%% ############## auxiliary functions for code readability ############################

% function that rotates vectors in the 2d plane by an angle theta
function [X,Y] = rotate2d(x,y,theta)
    X = x.*cos(theta) - y.*sin(theta);
    Y = x.*sin(theta) + y.*cos(theta);
end

% objective function to fit river Height H0 and river slope dH to a multilooked version of the
% radargram
function F = ObjFunc2(x,river_points,CS1b,river,N_segments,equality_constraint)
%     H0 = x(1);
%     dH_lat = x(2);
%     dH_lon = x(3);
    if ~isempty(equality_constraint)
        x(equality_constraint(1)) = x(equality_constraint(2));
    end
    mean_lat = mean(river.Lat,'omitnan');
    mean_lon = mean(river.Lon,'omitnan');
    rg_buffer = 7.5;
    %rg_buffer = 10;
    % initialize the range resolution and refbin
    dr = 0.18974206202531646;
    refBin = 257;
    F = 0;
    SARdata = CS1b.SAR.data20;
    normF = 0;
    for i = 1:20:numel(river_points)
        H_wgs84 = zeros(size(river_points(i).Lat));
        for seg_idx = 1:N_segments
            H_wgs84 = H_wgs84 + (seg_idx==river_points(i).segmentID).*(x(3*(seg_idx-1)+1) + x(3*(seg_idx-1)+2)*(river_points(i).Lat-mean_lat)  + x(3*(seg_idx-1)+3)*(river_points(i).Lon-mean_lon));
        end
        % now we have to substract the elevation that is due to the cross-track
        % distance over a curved surface, and can then use this height to
        % calculate the range in the radargram
        yr = H_wgs84 - river_points(i).yr0;
        Hs = CS1b.GEO.H(i);
        xr = river_points(i).xr;
        R = Hs + xr.^2/(2*Hs) - CS1b.MEA.tracker_range(i) - yr;
        rg = R/dr+refBin;
        rg(isnan(rg))=[];

        for k = 1:numel(rg)/2
            if rg(2*k-1)>rg(2*k)
               rgs =  round(rg(2*k)-rg_buffer:rg(2*k-1)+rg_buffer);
            else
                rgs =  round(rg(2*k-1)-rg_buffer:rg(2*k)+rg_buffer);
            end
            rgs(rgs>300)=[]; % for now
            rgs(rgs<1)=[];

            F = F + sum(log10(SARdata(rgs,i)));
            normF = normF + numel(rgs);
        end
    end
    F=-F/normF;
end