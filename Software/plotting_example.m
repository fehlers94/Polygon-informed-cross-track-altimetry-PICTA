filename = '/home/fehelers/PhD Delft/Projects/hydrology river/Results/L2_Garonne_cor/S6A_P4_1A_HR______20210527T182205_20210527T191822_20220509T142357_3377_020_070_035_EUM__REP_NT_F06.SEN6_TUDelft.nc'
water_level = ncread(filename,"H_wgs84_avg_cor");
lon = ncread(filename,"Lon_avg");
lat = ncread(filename,"Lat_avg");
lon_sat = ncread(filename,"Lon_sat");
lat_sat = ncread(filename,"Lat_sat");
mask = ncread(filename,"mask_valid");
mask = (mask==1);

fig = figure('units','inch','position',[0,0,8,8],'visible','on');
set(gcf,'color','w');

% ############### plot the map with elevations of the river segments

geoplot(lat_sat,lon_sat,'LineWidth',2,'LineStyle','--','Color','k'); hold on
geoscatter(lat(:),lon(:),'k.');hold on
geoscatter(lat(mask),lon(mask),[],water_level(mask),'filled');hold on

geolimits([44.3 44.6],[0.05 0.16])
caxis([50,70])
geotickformat('dd')
%geobasemap colorterrain
%geobasemap none
colorbar()
colormap(rgb)
legend('satellite track','river points','river height')
