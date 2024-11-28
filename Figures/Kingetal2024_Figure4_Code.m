%% Figure 4 
%Code to plot Figure 4 in King et al., 2024. Inputs are the Sea Level
%Histories from the reference earth model outputted by the Sea Level Model.

%% Create GaussGrid lat/long
xGrid = 0:0.25:360;
yGrid = -90:0.25:90;
[xq,yq] = meshgrid(xGrid, yGrid);
maxdeg = 256; 
N = maxdeg; 
[x,w] = GaussQuad(N); % Gauss Quadrature points (x) and weights (w)
x_GL = acos(x)*180/pi - 90; % cos(x_GL) are Quadrature points for P(cos(theta))
lon_GL = linspace(0,360,2*N+1);
lon_GL = lon_GL(1:end-1);

%puts it on a meshgrid
[lon_out,lat_out] = meshgrid(lon_GL,x_GL);

%% Aurora Basin

%Input the SL grid matrix both relative to present and normalized. The
%matrices are the the difference (i.e., sea-level change) between the KM3 
%interglacial and M2 glacial.
RelativeSLDataAuroraBasin = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_AuroraBasin_diff.txt');
NormalizedSLMapAuroraBasin =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_AuroraBasin_diff_normal.txt');

minz = 0; %min z value you want to show
maxz = 2; %max z value you want to show
numIntervals = 20; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('div', 'RdBu',numIntervals); %generating a color table

figure
set(gcf,'color','w') %background color
m_proj('robinson','lon',[0 360],'lat',[-90 90]); %projection
m_contourf(lon_out,lat_out,NormalizedSLMapAuroraBasin(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor',CT(1,:),'FontSize',14);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized sea level amplitude (m)';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;


%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%Plot locations on global sea level grid
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6);

% Virgina, USA
m_plot((360-75.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Enewetak Atoll, Pacific Ocean
m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

%Calcuate the sea level change ocurring at each of the three locations
%%%New Zealand
lat = -39.6964; lon = 175.5241;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
ABNewZealandAmplitude = squeeze(RelativeSLDataAuroraBasin(ilat,ilon,:));

%%%Virginia
lat = 37.321473; lon = 284.0248;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
ABVirginiaAmplitude = squeeze(RelativeSLDataAuroraBasin(ilat,ilon,:));

%%%Enewetak Atoll
lat = 11.50; lon = 162.33;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
ABEnewetakAtollAmplitude = squeeze(RelativeSLDataAuroraBasin(ilat,ilon,:));

%% East Antarctica
RelativeSLDataEastAntarctica = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_EastAntarctica_diff.txt');
NormalizedSLMapEastAntarctica =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_EastAntarctica_diff_normal.txt');

minz = 0; %min z value you want to show
maxz = 2; %max z value you want to show
numIntervals = 20; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('div', 'RdBu',numIntervals); %generating a color table
%CM=cmocean('amp',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,NormalizedSLMapEastAntarctica(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor',CT(1,:),'FontSize',14);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized sea level amplitude (m)';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Virgina, USA
m_plot((360-75.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Enewetak Atoll, Pacific Ocean
m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

%%%New Zealand
lat = -39.6964; lon = 175.5241;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
EAISNewZealandAmplitude = squeeze(RelativeSLDataEastAntarctica(ilat,ilon,:));

%%%Virginia
lat = 37.321473; lon = 284.0248;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
EAISVirginiaAmplitude = squeeze(RelativeSLDataEastAntarctica(ilat,ilon,:));

%%%Enewetak Atoll
lat = 11.50; lon = 162.33;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
EAISEnewetakAtollAmplitude = squeeze(RelativeSLDataEastAntarctica(ilat,ilon,:));

%% Eurasia
RelativeSLDataEurasia = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_Eurasia_diff.txt');
NormalizedSLMapEurasia =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_Eurasia_diff_normal.txt');


minz = 0; %min z value you want to show
maxz = 2; %max z value you want to show
numIntervals = 20; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('div', 'RdBu',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,NormalizedSLMapEurasia(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor',CT(1,:),'FontSize',14);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized sea level amplitude (m)';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;


%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Virgina, USA
m_plot((360-71.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Enewetak Atoll, Pacific Ocean
m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

%%%New Zealand
lat = -39.6964; lon = 175.5241;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
EISNewZealandAmplitude = squeeze(RelativeSLDataEurasia(ilat,ilon,:));

%%%Virginia
lat = 37.321473; lon = 284.0248;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
EISVirginiaAmplitude = squeeze(RelativeSLDataEurasia(ilat,ilon,:));

%%%Enewetak Atoll
lat = 11.50; lon = 162.33;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
EISEnewetakAtollAmplitude = squeeze(RelativeSLDataEurasia(ilat,ilon,:));

%% Greenland
RelativeSLDataGreenland = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_Greenland_diff.txt');
NormalizedSLMapGreenland =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_Greenland_diff_normal.txt');

minz = 0; %min z value you want to show
maxz = 2; %max z value you want to show
numIntervals = 20; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('div', 'RdBu',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,NormalizedSLMapGreenland(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor',CT(1,:),'FontSize',14);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized sea level amplitude (m)';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Virgina, USA
m_plot((360-75.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Enewetak Atoll, Pacific Ocean
m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

%%%New Zealand
lat = -39.6964; lon = 175.5241;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
GrISNewZealandAmplitude = squeeze(RelativeSLDataGreenland(ilat,ilon,:));

%%%Virginia
lat = 37.321473; lon = 284.0248;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
GrISVirginiaAmplitude = squeeze(RelativeSLDataGreenland(ilat,ilon,:));

%%%Enewetak Atoll
lat = 11.50; lon = 162.33;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
GrISEnewetakAtollAmplitude = squeeze(RelativeSLDataGreenland(ilat,ilon,:));

%%  North  America (35 m)
RelativeSLDataNorthAmerica = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_NorthAmerica35m_diff.txt');
NormalizedSLMapNorthAmerica =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_NorthAmerica35m_diff_normal.txt');

minz = 0; %min z value you want to show
maxz = 2; %max z value you want to show
numIntervals = 20; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('div', 'RdBu',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,NormalizedSLMapNorthAmerica(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor',CT(1,:),'FontSize',14);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized sea level amplitude (m)';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Virgina, USA
m_plot((360-75.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Enewetak Atoll, Pacific Ocean
m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

%%%New Zealand
lat = -39.6964; lon = 175.5241;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
NAISNewZealandAmplitude = squeeze(RelativeSLDataNorthAmerica(ilat,ilon,:));

%%%Virginia
lat = 37.321473; lon = 284.0248;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
NAISVirginiaAmplitude = squeeze(RelativeSLDataNorthAmerica(ilat,ilon,:));

%%%Enewetak Atoll
lat = 11.50; lon = 162.33;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
NAISEnewetakAtollAmplitude = squeeze(RelativeSLDataNorthAmerica(ilat,ilon,:));
%% North  America (8 m)
RelativeSLDataNorthAmerica = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_NorthAmerica8m_diff.txt');
NormalizedSLMapNorthAmerica =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_NorthAmerica8m_diff_normal.txt');

minz = 0; %min z value you want to show
maxz = 2; %max z value you want to show
numIntervals = 20; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('div', 'RdBu',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,NormalizedSLMapNorthAmerica(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor',CT(1,:),'FontSize',14);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized sea level amplitude (m)';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Virgina, USA
m_plot((360-75.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Enewetak Atoll, Pacific Ocean
m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

%%%New Zealand
lat = -39.6964; lon = 175.5241;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
NAISNewZealandAmplitude = squeeze(RelativeSLDataNorthAmerica(ilat,ilon,:));

%%%Virginia
lat = 37.321473; lon = 284.0248;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
NAISVirginiaAmplitude = squeeze(RelativeSLDataNorthAmerica(ilat,ilon,:));

%%%Enewetak Atoll
lat = 11.50; lon = 162.33;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
NAISEnewetakAtollAmplitude = squeeze(RelativeSLDataNorthAmerica(ilat,ilon,:));

%% Prydz Bay
RelativeSLDataPrydzBay = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_PrydzBay_diff.txt');
NormalizedSLMapPrydzBay =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_PrydzBay_diff_normal.txt');

minz = 0; %min z value you want to show
maxz = 2; %max z value you want to show
numIntervals = 20; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('div', 'RdBu',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,NormalizedSLMapPrydzBay(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor',CT(1,:),'FontSize',14);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized sea level amplitude (m)';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Virgina, USA
m_plot((360-75.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Enewetak Atoll, Pacific Ocean
m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

%%%New Zealand
lat = -39.6964; lon = 175.5241;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
PBNewZealandAmplitude = squeeze(RelativeSLDataPrydzBay(ilat,ilon,:));

%%%Virginia
lat = 37.321473; lon = 284.0248;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
PBVirginiaAmplitude = squeeze(RelativeSLDataPrydzBay(ilat,ilon,:));

%%%Enewetak Atoll
lat = 11.50; lon = 162.33;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
PBEnewetakAtollAmplitude = squeeze(RelativeSLDataPrydzBay(ilat,ilon,:));

%% West Antarctica
RelativeSLDataWestAntarctica = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_WestAntarctica_diff.txt');
NormalizedSLMapWestAntarctica =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_WestAntarctica_diff_normal.txt');

minz = 0; %min z value you want to show
maxz = 2; %max z value you want to show
numIntervals = 20; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('div', 'RdBu',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,NormalizedSLMapWestAntarctica(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor',CT(1,:),'FontSize',14);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized sea level amplitude (m)';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Virgina, USA
m_plot((360-75.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Enewetak Atoll, Pacific Ocean
m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

%%%New Zealand
lat = -39.6964; lon = 175.5241;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
WAISNewZealandAmplitude = squeeze(RelativeSLDataWestAntarctica(ilat,ilon,:));

%%%Virginia
lat = 37.321473; lon = 284.0248;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
WAISVirginiaAmplitude = squeeze(RelativeSLDataWestAntarctica(ilat,ilon,:));

%%%Enewetak Atoll
lat = 11.50; lon = 162.33;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
WAISEnewetakAtollAmplitude = squeeze(RelativeSLDataWestAntarctica(ilat,ilon,:));

%% Wilkes Basin
RelativeSLDataWilkesBasin = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_WilkesBasin_diff.txt');
NormalizedSLMapWilkesBasin =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_WilkesBasin_diff_normal.txt');

minz = 0; %min z value you want to show
maxz = 2; %max z value you want to show
numIntervals = 20; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('div', 'RdBu',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,NormalizedSLMapWilkesBasin(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor',CT(1,:),'FontSize',14);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized sea level amplitude (m)';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Virgina, USA
m_plot((360-75.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

% Enewetak Atoll, Pacific Ocean
m = m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on

%%%New Zealand
lat = -39.6964; lon = 175.5241;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
WBNewZealandAmplitude = squeeze(RelativeSLDataWilkesBasin(ilat,ilon,:));

%%%Virginia
lat = 37.321473; lon = 284.0248;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
WBVirginiaAmplitude = squeeze(RelativeSLDataWilkesBasin(ilat,ilon,:));

%%%Enewetak Atoll
lat = 11.50; lon = 162.33;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
WBEnewetakAtollAmplitude = squeeze(RelativeSLDataWilkesBasin(ilat,ilon,:));