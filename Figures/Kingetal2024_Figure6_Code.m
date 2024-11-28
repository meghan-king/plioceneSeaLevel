%% Figure 6
%Code to plot Figure 6 in King et al., 2024. Inputs are the Normalized Sea Level
%Histories from the reference earth model outputted by the Sea Level Model.
%The matrices are the the difference (i.e., sea-level change) between the KM3 
%interglacial and M2 glacial.

%% Input Normalized Sea Level Grids
NormalizedSLMapAuroraBasin =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_AuroraBasin_diff_normal.txt');
NormalizedSLMapEastAntarctica =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_EastAntarctica_diff_normal.txt');
NormalizedSLMapEurasia =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_Eurasia_diff_normal.txt');
NormalizedSLMapGreenland =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_Greenland_diff_normal.txt');
NormalizedSLMapNorthAmerica =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_NorthAmerica35m_diff_normal.txt');
NormalizedSLMapNorthAmerica8m =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_NorthAmerica8m_diff_normal.txt');
NormalizedSLMapPrydzBay =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_PrydzBay_diff_normal.txt');
NormalizedSLMapWestAntarctica =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_WestAntarctica_diff_normal.txt');
NormalizedSLMapWilkesBasin =  readmatrix('.../SEAGL_GRIDS/seagl_96Cp55_WilkesBasin_diff_normal.txt');

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

%% Fig 6a: maximum percentage disrepancy from GMSL based on a contribution from East Antarctica, Greenland, and West Antarctica 

%Calcuate areas that are within 5% of GMSL
EastAntarctica5 = NormalizedSLMapEastAntarctica;
EastAntarctica5(EastAntarctica5>0.95 & EastAntarctica5<1.05)= .03;
EastAntarctica5(EastAntarctica5 ~=.03)= NaN;

Greenland5 = NormalizedSLMapGreenland;
Greenland5(Greenland5>0.95 & Greenland5<1.05)= .01;
Greenland5(Greenland5 ~=.01)= NaN;

WestAntarctica5 = NormalizedSLMapWestAntarctica;
WestAntarctica5(WestAntarctica5>0.95 & WestAntarctica5<1.05)= .01;
WestAntarctica5(WestAntarctica5 ~=.01)= NaN;

NormalizedSLMapCombinationIce5 = (EastAntarctica5 + Greenland5 + WestAntarctica5);
NormalizedSLMapCombinationIce5(isnan(NormalizedSLMapCombinationIce5)) = 0;

%Calcuate areas that are within 10% of GMSL
EastAntarctica10 = NormalizedSLMapEastAntarctica;
EastAntarctica10(EastAntarctica10>0.9 & EastAntarctica10<1.1)= .06;
EastAntarctica10(EastAntarctica10 ~=.06)= NaN;

Greenland10 = NormalizedSLMapGreenland;
Greenland10(Greenland10>0.9 & Greenland10<1.1)= .02;
Greenland10(Greenland10 ~=.02)= NaN;

WestAntarctica10 = NormalizedSLMapWestAntarctica;
WestAntarctica10(WestAntarctica10>0.9 & WestAntarctica10<1.1)= .02;
WestAntarctica10(WestAntarctica10 ~=.02)= NaN;

NormalizedSLMapCombinationIce10 = (EastAntarctica10 + Greenland10 + WestAntarctica10);
NormalizedSLMapCombinationIce10(isnan(NormalizedSLMapCombinationIce10)) = 0;

%Calcuate areas that are within 15% of GMSL
EastAntarctica15 = NormalizedSLMapEastAntarctica;
EastAntarctica15(EastAntarctica15>0.85 & EastAntarctica15<1.15)= .09;
EastAntarctica15(EastAntarctica15 ~=.09)= NaN;

Greenland15 = NormalizedSLMapGreenland;
Greenland15(Greenland15>0.85 & Greenland15<1.15)= .03;
Greenland15(Greenland15 ~=.03)= NaN;

WestAntarctica15 = NormalizedSLMapWestAntarctica;
WestAntarctica15(WestAntarctica15>0.85 & WestAntarctica15<1.15)= .03;
WestAntarctica15(WestAntarctica15 ~=.03)= NaN;

NormalizedSLMapCombinationIce15 = (EastAntarctica15 + Greenland15 + WestAntarctica15);
NormalizedSLMapCombinationIce15(isnan(NormalizedSLMapCombinationIce15)) = 0;

%Calcuate areas that are within 20% of GMSL
EastAntarctica20 = NormalizedSLMapEastAntarctica;
EastAntarctica20(EastAntarctica20>0.80 & EastAntarctica20<1.20)= .12;
EastAntarctica20(EastAntarctica20 ~=.12)= NaN;

Greenland20 = NormalizedSLMapGreenland;
Greenland20(Greenland20>0.80 & Greenland20<1.20)= .04;
Greenland20(Greenland20 ~=.04)= NaN;

WestAntarctica20 = NormalizedSLMapWestAntarctica;
WestAntarctica20(WestAntarctica20>0.80 & WestAntarctica20<1.20)= .04;
WestAntarctica20(WestAntarctica20 ~=.04)= NaN;

NormalizedSLMapCombinationIce20 = (EastAntarctica20 + Greenland20 + WestAntarctica20);
NormalizedSLMapCombinationIce20(isnan(NormalizedSLMapCombinationIce20)) = 0;

%Calcuate areas that are 25% of GMSL or greater
EastAntarctica25 = NormalizedSLMapEastAntarctica;
EastAntarctica25(EastAntarctica25>0.75 & EastAntarctica25<1.25)= .15;
EastAntarctica25(EastAntarctica25 ~=.15)= NaN;

Greenland25 = NormalizedSLMapGreenland;
Greenland25(Greenland25>0.75 & Greenland25<1.25)= .05;
Greenland25(Greenland25 ~=.05)= NaN;

WestAntarctica25 = NormalizedSLMapWestAntarctica;
WestAntarctica25(WestAntarctica25>0.75 & WestAntarctica25<1.25)= .05;
WestAntarctica25(WestAntarctica25 ~=.05)= NaN;

NormalizedSLMapCombinationIce25 = (EastAntarctica25 + Greenland25 + WestAntarctica25);
NormalizedSLMapCombinationIce25(isnan(NormalizedSLMapCombinationIce25)) = 0;

%Final Map
%Combine these areas onto a global grid and plot each percentage as a
%different shade of purple
FinalMap = NormalizedSLMapCombinationIce25+NormalizedSLMapCombinationIce20;
FinalMap(FinalMap==0.45) = 0.20;

FinalMap = FinalMap+ NormalizedSLMapCombinationIce15;
FinalMap(FinalMap==0.35) = 0.15;

FinalMap = FinalMap+ (NormalizedSLMapCombinationIce10*2);
FinalMap(FinalMap==0.35) = 0.10;

FinalMap = FinalMap+ (NormalizedSLMapCombinationIce5*5);
FinalMap(FinalMap==0.35) = 0.05;


FinalMap(FinalMap ==0)= NaN;
FinalMap(FinalMap ==0.25)= 0.24;
FinalMap(FinalMap ==0.20)= 0.19;
FinalMap(FinalMap ==0.15)= 0.14;
FinalMap(FinalMap ==0.10)= 0.09;
FinalMap(FinalMap ==0.05)= 0.04;
 

minz = 0; %min z value you want to show
maxz = 0.25; %max z value you want to show
numIntervals = 5; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('seq', 'Purples',numIntervals); %generating a color table
CT = flip(CT);
%CM=cmocean('amp',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,FinalMap(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',16); %draw in a grid
m_grid('backgroundcolor','w','FontSize',16);
% m_coast(); %add in continents
% m_coast('patch',[.9 .9 .9],'EdgeColor',[.75 .75 .75]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized GMSL';
cb.Label.FontSize = 16;
% cb.TickLabels.FontSize = 14;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); %Add point to map (10-30m;20?;Naish & Wilson, 2009; Grant et al., 2019)
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); %Add point to map (Naish & Wilson, 2009; Grant et al., 2019)

% Virgina, USA
m_plot((360-75.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); %Add point to map (15-20m;Krantz, 1991)

% Enewetak Atoll, Pacific Ocean
m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); %Add point to map (Wardlaw & Quinn, 1991)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 6b: maximum percentage disrepancy from GMSL based on a contribution from North America, East Antarctica, Greenland, West Antarctica and Eurasia

%5
EastAntarctica5 = NormalizedSLMapEastAntarctica;
EastAntarctica5(EastAntarctica5>0.95 & EastAntarctica5<1.05)= .01;
EastAntarctica5(EastAntarctica5 ~=.01)= NaN;

Eurasia5 = NormalizedSLMapEurasia;
Eurasia5(Eurasia5>0.95 & Eurasia5<1.05)= .01;
Eurasia5(Eurasia5 ~=.01)= NaN;

Greenland5 = NormalizedSLMapGreenland;
Greenland5(Greenland5>0.95 & Greenland5<1.05)= .01;
Greenland5(Greenland5 ~=.01)= NaN;

NorthAmerica5 = NormalizedSLMapNorthAmerica;
NorthAmerica5(NorthAmerica5>0.95 & NorthAmerica5<1.05)= .01;
NorthAmerica5(NorthAmerica5 ~=.01)= NaN;

WestAntarctica5 = NormalizedSLMapWestAntarctica;
WestAntarctica5(WestAntarctica5>0.95 & WestAntarctica5<1.05)= .01;
WestAntarctica5(WestAntarctica5 ~=.01)= NaN;

NormalizedSLMapCombinationIce5 = (EastAntarctica5 + Eurasia5 + Greenland5 + NorthAmerica5 + WestAntarctica5);
NormalizedSLMapCombinationIce5(isnan(NormalizedSLMapCombinationIce5)) = 0;

%10
EastAntarctica10 = NormalizedSLMapEastAntarctica;
EastAntarctica10(EastAntarctica10>0.9 & EastAntarctica10<1.1)= .02;
EastAntarctica10(EastAntarctica10 ~=.02)= NaN;

Eurasia10 = NormalizedSLMapEurasia;
Eurasia10(Eurasia10>0.9 & Eurasia10<1.1)= .02;
Eurasia10(Eurasia10 ~=.02)= NaN;

Greenland10 = NormalizedSLMapGreenland;
Greenland10(Greenland10>0.9 & Greenland10<1.1)= .02;
Greenland10(Greenland10 ~=.02)= NaN;

NorthAmerica10 = NormalizedSLMapNorthAmerica;
NorthAmerica10(NorthAmerica10>0.9 & NorthAmerica10<1.1)= .02;
NorthAmerica10(NorthAmerica10 ~=.02)= NaN;

WestAntarctica10 = NormalizedSLMapWestAntarctica;
WestAntarctica10(WestAntarctica10>0.9 & WestAntarctica10<1.1)= .02;
WestAntarctica10(WestAntarctica10 ~=.02)= NaN;

NormalizedSLMapCombinationIce10 = (EastAntarctica10 + Eurasia10 + Greenland10 + NorthAmerica10 + WestAntarctica10);
NormalizedSLMapCombinationIce10(isnan(NormalizedSLMapCombinationIce10)) = 0;

%15
EastAntarctica15 = NormalizedSLMapEastAntarctica;
EastAntarctica15(EastAntarctica15>0.85 & EastAntarctica15<1.15)= .03;
EastAntarctica15(EastAntarctica15 ~=.03)= NaN;

Eurasia15 = NormalizedSLMapEurasia;
Eurasia15(Eurasia15>0.85 & Eurasia15<1.15)= .03;
Eurasia15(Eurasia15 ~=.03)= NaN;

Greenland15 = NormalizedSLMapGreenland;
Greenland15(Greenland15>0.85 & Greenland15<1.15)= .03;
Greenland15(Greenland15 ~=.03)= NaN;

NorthAmerica15 = NormalizedSLMapNorthAmerica;
NorthAmerica15(NorthAmerica15>0.85 & NorthAmerica15<1.15)= .03;
NorthAmerica15(NorthAmerica15 ~=.03)= NaN;

WestAntarctica15 = NormalizedSLMapWestAntarctica;
WestAntarctica15(WestAntarctica15>0.85 & WestAntarctica15<1.15)= .03;
WestAntarctica15(WestAntarctica15 ~=.03)= NaN;

NormalizedSLMapCombinationIce15 = (EastAntarctica15 + Eurasia15 + Greenland15 + NorthAmerica15 + WestAntarctica15);
NormalizedSLMapCombinationIce15(isnan(NormalizedSLMapCombinationIce15)) = 0;

%20
EastAntarctica20 = NormalizedSLMapEastAntarctica;
EastAntarctica20(EastAntarctica20>0.80 & EastAntarctica20<1.20)= .04;
EastAntarctica20(EastAntarctica20 ~=.04)= NaN;

Eurasia20 = NormalizedSLMapEurasia;
Eurasia20(Eurasia20>0.80 & Eurasia20<1.20)= .04;
Eurasia20(Eurasia20 ~=.04)= NaN;

Greenland20 = NormalizedSLMapGreenland;
Greenland20(Greenland20>0.80 & Greenland20<1.20)= .04;
Greenland20(Greenland20 ~=.04)= NaN;

NorthAmerica20 = NormalizedSLMapNorthAmerica;
NorthAmerica20(NorthAmerica20>0.80 & NorthAmerica20<1.20)= .04;
NorthAmerica20(NorthAmerica20 ~=.04)= NaN;

WestAntarctica20 = NormalizedSLMapWestAntarctica;
WestAntarctica20(WestAntarctica20>0.80 & WestAntarctica20<1.20)= .04;
WestAntarctica20(WestAntarctica20 ~=.04)= NaN;

NormalizedSLMapCombinationIce20 = (EastAntarctica20 + Eurasia20 + Greenland20 + NorthAmerica20 + WestAntarctica20);
NormalizedSLMapCombinationIce20(isnan(NormalizedSLMapCombinationIce20)) = 0;

%25
EastAntarctica25 = NormalizedSLMapEastAntarctica;
EastAntarctica25(EastAntarctica25>0.75 & EastAntarctica25<1.25)= .05;
EastAntarctica25(EastAntarctica25 ~=.05)= NaN;

Eurasia25 = NormalizedSLMapEurasia;
Eurasia25(Eurasia25>0.75 & Eurasia25<1.25)= .05;
Eurasia25(Eurasia25 ~=.05)= NaN;

Greenland25 = NormalizedSLMapGreenland;
Greenland25(Greenland25>0.75 & Greenland25<1.25)= .05;
Greenland25(Greenland25 ~=.05)= NaN;

NorthAmerica25 = NormalizedSLMapNorthAmerica;
NorthAmerica25(NorthAmerica25>0.75 & NorthAmerica25<1.25)= .05;
NorthAmerica25(NorthAmerica25 ~=.05)= NaN;

WestAntarctica25 = NormalizedSLMapWestAntarctica;
WestAntarctica25(WestAntarctica25>0.75 & WestAntarctica25<1.25)= .05;
WestAntarctica25(WestAntarctica25 ~=.05)= NaN;

NormalizedSLMapCombinationIce25 = (EastAntarctica25 + Eurasia25 + Greenland25 + NorthAmerica25 + WestAntarctica25);
NormalizedSLMapCombinationIce25(isnan(NormalizedSLMapCombinationIce25)) = 0;

%Final Map
FinalMap = NormalizedSLMapCombinationIce25+NormalizedSLMapCombinationIce20;
FinalMap(FinalMap==0.45) = 0.20;

FinalMap = FinalMap+ NormalizedSLMapCombinationIce15;
FinalMap(FinalMap==0.35) = 0.15;

FinalMap = FinalMap+ (NormalizedSLMapCombinationIce10*2);
FinalMap(FinalMap==0.35) = 0.10;

FinalMap = FinalMap+ (NormalizedSLMapCombinationIce5*5);
FinalMap(FinalMap==0.35) = 0.05;


FinalMap(FinalMap ==0)= NaN;
FinalMap(FinalMap ==0.25)= 0.24;
FinalMap(FinalMap ==0.20)= 0.19;
FinalMap(FinalMap ==0.15)= 0.14;
FinalMap(FinalMap ==0.10)= 0.09;
FinalMap(FinalMap ==0.05)= 0.04;
 

minz = 0.0; %min z value you want to show
maxz = 0.25; %max z value you want to show
numIntervals = 5; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('seq', 'Purples',numIntervals); %generating a color table
CT = flip(CT);
%CM=cmocean('amp',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_contourf(lon_out,lat_out,FinalMap(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',16); %draw in a grid
m_grid('backgroundcolor','w','FontSize',16);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.75 .75 .75]); %add in continents
cb = colorbar;
cb.Label.String = 'Normalized GMSL';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%%%%%%%%%%%%%%% Physical Records of Pliocene Sea Level %%%%%%%%%%%%%%%%%%%
%%% Whanganui Baisn, New Zealand
m_plot(175.5241,-39.6964,'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); %Add point to map (10-30m;20?;Naish & Wilson, 2009; Grant et al., 2019)
m_plot((175.807002),-39.812240, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); %Add point to map (Naish & Wilson, 2009; Grant et al., 2019)

% Virgina, USA
m_plot((360-75.98),37.32, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); %Add point to map (15-20m;Krantz, 1991)

% Enewetak Atoll, Pacific Ocean
m_plot((162.33),11.50, 'wo', 'MarkerFaceColor', 'k', 'MarkerSize',6); %Add point to map (Wardlaw & Quinn, 1991)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
