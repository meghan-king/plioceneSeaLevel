%% Figure 3 
%Code to plot Figure 3 in King et al., 2024. Inputs are the Ice Histories
%outputted by the Ice Model. 

%% AuroraBasin
data = fopen('.../IceModels/AuroraBasin/ice316'); %Time slice ice316 is the M2 glaciation
M2AuroraBasin = fread(data, 'float32'); 
data = fopen('.../IceModels/AuroraBasin/ice456'); %Time slice ice456 is the KM3 glaciation
KM3AuroraBasin = fread(data, 'float32');

if length(M2AuroraBasin)== 131072 %Ice model output comes out like this sometimes
        M2AuroraBasin = reshape(M2AuroraBasin,512,256)'; %Reshape Guass Grid into 256x512
        KM3AuroraBasin = reshape(KM3AuroraBasin,512,256)'; %Reshape Guass Grid into 256x512
end

%% EastAntarctica
data = fopen('.../IceModels/EastAntarctica/ice316'); %131,346,316
M2EastAntarctica = fread(data, 'float32');
data = fopen('.../IceModels/EastAntarctica/ice456'); %131,346,316
KM3EastAntarctica = fread(data, 'float32');

if length(M2EastAntarctica)== 131072 %Ice model output comes out like this sometimes
        M2EastAntarctica = reshape(M2EastAntarctica,512,256)';
        KM3EastAntarctica = reshape(KM3EastAntarctica,512,256)';
end

%% Eurasia
data = fopen('.../IceModels/Eurasia/ice316'); %131,346,316
M2Eurasia = fread(data, 'float32');
data = fopen('.../IceModels/Eurasia/ice456'); %131,346,316
KM3Eurasia = fread(data, 'float32');

if length(M2Eurasia)== 131072 %Ice model output comes out like this sometimes
        M2Eurasia = reshape(M2Eurasia,512,256)';
        KM3Eurasia = reshape(KM3Eurasia,512,256)';
end

%% Greenland
data = fopen('/Volumes/Extreme SSD/IceModels/Greenland/ice316'); %131,346,316
M2Greenland = fread(data, 'float32');
data = fopen('/Volumes/Extreme SSD/IceModels/Greenland/ice456'); %131,346,316
KM3Greenland = fread(data, 'float32');

if length(M2Greenland)== 131072 %Ice model output comes out like this sometimes
        M2Greenland = reshape(M2Greenland,512,256)';
        KM3Greenland = reshape(KM3Greenland,512,256)';
end

%% NorthAmerica
data = fopen('/Volumes/Extreme SSD/IceModels/NorthAmerica/ice316'); %131,346,316
M2NorthAmerica = fread(data, 'float32');
data = fopen('/Volumes/Extreme SSD/IceModels/NorthAmerica/ice456'); %131,346,316
KM3NorthAmerica = fread(data, 'float32');

if length(M2NorthAmerica)== 131072 %Ice model output comes out like this sometimes
        M2NorthAmerica = reshape(M2NorthAmerica,512,256)';
        KM3NorthAmerica = reshape(KM3NorthAmerica,512,256)';
end

%% PrydzBay
data = fopen('/Volumes/Extreme SSD/IceModels/PrydzBay/ice316'); %131,346,316
M2PrydzBay = fread(data, 'float32');
data = fopen('/Volumes/Extreme SSD/IceModels/PrydzBay/ice456'); %131,346,316
KM3PrydzBay = fread(data, 'float32');

if length(M2PrydzBay)== 131072 %Ice model output comes out like this sometimes
        M2PrydzBay = reshape(M2PrydzBay,512,256)';
        KM3PrydzBay = reshape(KM3PrydzBay,512,256)';
end

%% WestAntarctica
data = fopen('/Volumes/Extreme SSD/PhD Research/Chapter 2/Pliocene Data/Sea Level Model Data/ExtraIceModels/WestAntarctica/ice316'); %131,346,316
M2WestAntarctica = fread(data, 'float32');
data = fopen('/Volumes/Extreme SSD/PhD Research/Chapter 2/Pliocene Data/Sea Level Model Data/ExtraIceModels/WestAntarctica/ice456'); %131,346,316
KM3WestAntarctica = fread(data, 'float32');

if length(M2WestAntarctica)== 131072 %Ice model output comes out like this sometimes
        M2WestAntarctica = reshape(M2WestAntarctica,512,256)';
        KM3WestAntarctica = reshape(KM3WestAntarctica,512,256)';
end

%% WilkesBasin
data = fopen('/Volumes/Extreme SSD/IceModels/WilkesBasin/ice316'); %131,346,316
M2WilkesBasin = fread(data, 'float32');
data = fopen('/Volumes/Extreme SSD/IceModels/WilkesBasin/ice456'); %131,346,316
KM3WilkesBasin = fread(data, 'float32');

if length(M2WilkesBasin)== 131072 %Ice model output comes out like this sometimes
        M2WilkesBasin = reshape(M2WilkesBasin,512,256)';
        KM3WilkesBasin = reshape(KM3WilkesBasin,512,256)';
end

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

%Puts it on a meshgri
[lon_out,lat_out] = meshgrid(lon_GL,x_GL);

%% Create M2 glacial global ice plot

%Combine indvidual M2 ice files
M2CombinedIce =  M2EastAntarctica + M2Eurasia + M2Greenland + M2NorthAmerica  + M2WestAntarctica ;
M2CombinedIce(M2CombinedIce==0)=NaN;

minz = 0; %min z value you want to show
maxz = 1000; %max z value you want to show
numIntervals = 10; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('seq', 'Blues',numIntervals); %generating a color table
%CM=cmocean('amp',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
hold on;
m_contourf(lon_out,lat_out,M2CombinedIce(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor','none','FontSize',14);

cb = colorbar;
cb.Label.String = 'Ice Thickness (m)';
cb.Label.FontSize = 16;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;

%% Create KM3 Interglacial Global Ice Plot

%Combine indvidual KM3 ice files
KM3CombinedIce =  KM3EastAntarctica + KM3Eurasia + KM3Greenland + KM3NorthAmerica + KM3WestAntarctica ;
KM3CombinedIce(KM3CombinedIce==0)=NaN;

minz = 0; %min z value you want to show
maxz = 1000; %max z value you want to show
numIntervals = 10; %This is the number of color slices you'll have
dL = (maxz-minz)/numIntervals; %this is the interval size
L=[minz:dL:maxz]; %these are the contour levels
CT=cbrewer('seq', 'Blues',numIntervals); %generating a color table
%CM=cmocean('amp',numIntervals); %generating a color table

figure
set(gcf,'color','w')
m_proj('robinson','lon',[0 360],'lat',[-90 90]);
m_coast(); %add in continents
m_coast('patch',[.9 .9 .9],'EdgeColor',[.6 .6 .6]); %add in continents
hold on;
m_contourf(lon_out,lat_out,KM3CombinedIce(:,:),L,'linecolor','none'); %contour plot 
hold on;
colormap(CT) %set colormap of this figure to the one made by cbrewer
caxis([minz maxz])
colorbar('Ticks',L, 'TickLabels', L,'FontSize',14); %colorbar with ticks and tick labels at contour levels
m_grid('tickdir','in','xticklabels',[],'yticklabels',[],'FontSize',14); %draw in a grid
m_grid('backgroundcolor','none','FontSize',14);

cb = colorbar;
cb.Label.String = 'Ice Thickness (m)';
cb.Label.FontSize = 16;
% cb.TickLabels.FontSize = 14;
a = gca;
a.LineWidth = 2;
a.Position(3) = 0.75;
% annotation('textbox', [.89, .01, .1, .1], 'String', "Underestimation",'EdgeColor','w','FontSize',14)
% annotation('textbox', [.89, .85, .1, .1], 'String', "Overestimation",'EdgeColor','w','FontSize',14)
% title("Normalized Sea Level Change (M2 to KM3) - East Antarctica",'FontSize',18)