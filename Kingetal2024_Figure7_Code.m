%% Figure 7
%Code to plot Figure 7 in King et al., 2024. Inputs are the Relative Sea Level
%Histories from the reference earth model outputted by the Sea Level Model.
%The matrices are the the difference (i.e., sea-level change) between the KM3 
%interglacial and M2 glacial.

%% Input
RelativeSLDataEastAntarctica = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_EastAntarctica_diff.txt');
RelativeSLDataEurasia = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_Eurasia_diff.txt');
RelativeSLDataGreenland = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_Greenland_diff.txt');
RelativeSLDataNorthAmerica = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_NorthAmerica35m_diff.txt');
RelativeSLDataNorthAmerica8m = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_NorthAmerica8m_diff.txt');
RelativeSLDataWestAntarctica = readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_WestAntarctica_diff.txt');

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
%% Calculate percentages of ice sheets to be combined
%Percentage of ice sheet to be factored into combination scenarios. These
%values can be changed based on the amount of each ice sheet to be
%contributed to the global sea level grid.
NAIS_value = 0.05;
GrIS_value = 1;
EIS_value = 0.5;
EAIS_value = 0.345;
WAIS_value = 1;

%GMSLP values for each ice sheet based on the reference earth model. Do not
%change.
NAIS_GMSLP = 33.12;
GrIS_GMSLP = 6.67;
EIS_GMSLP = 4.23;
EAIS_GMSLP = 11.08;
WAIS_GMSLP = 2.74;

%Scale sea level grids for each ice sheet based on the inputted _value, and then combine them into one plot. 
SL_CombinationPlot =  ((RelativeSLDataNorthAmerica)*NAIS_value) + ((RelativeSLDataGreenland)*GrIS_value) + ((RelativeSLDataEurasia)*EIS_value) + ((RelativeSLDataEastAntarctica)*EAIS_value) + ((RelativeSLDataWestAntarctica)*WAIS_value);

%Calcuate the global sea-level change. 
GMSLAmplitudeP = ((NAIS_GMSLP*NAIS_value) + (GrIS_GMSLP*GrIS_value) + (EIS_GMSLP*EIS_value) + (EAIS_GMSLP*EAIS_value) + (WAIS_GMSLP*WAIS_value));

%% Create Combination Plots
minz = -50; %min z value you want to show
maxz = 50; %max z value you want to show
numIntervals = 20; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('div', 'RdBu',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,SL_CombinationPlot(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',16); %draw in a grid
m_grid('backgroundcolor',CT(1,:),'FontSize',16);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor','none'); %add in continents
cb = colorbar;
cb.Label.String = 'sea level amplitude (m)';
cb.Label.FontSize = 16;
% cb.TickLabels.FontSize = 14;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calcuate the sea-level change ocurring at New Zealand
lat = -39.6964; lon = 175.5241;
lat_256 = lat_out(:,1); lon_256 = lon_out(1,:);
[~, ilat] = min(abs(lat_256 - lat));
[~, ilon] = min(abs(lon_256 - lon));
aNewZealandAmplitude = squeeze(SL_CombinationPlot(ilat,ilon,:)); 