%% Ice Model code for King et al., 2024
%This code creates individual ice sheet histories for the time period of
%3.61 - 2.95 Ma using snapshots of ice geometries from Berends et al., 2019 and Dowsett 
%The inputs are a d18O time series .txt file and ca MATLAB workspace .mat 
%file that contains snapshots of ice geometries. The outputs are binary files in 
%a 512x256 Gaussian Grrid containing an individual ice sheet geometry every
%1000 years (for 661,000 years) based on the d18O time series. The code was written in MATLAB 
%R2019a and may not work with newer versions of the software.

%% Part 1: Create normalized Lisiecki & Raymo time series

Input_data = importdata('Kingetal2024_IceModel_Inputd18O.txt'); %Input raw Lisiecki & Raymo d18O data

Time = flip(Input_data(:,1)); %Create Time vector from input
d18O = flip(Input_data(:,2)); %Create d18O vector from input

neg_d18O = (d18O)*-1; %Get the negative of the d18O curve
min_d18O = min(neg_d18O); %Create min d18O vector
max_d18O = max(neg_d18O); %Create max d18) vector
Scalar_factor_initial = zeros(length(Time),1); %Create scalar_factor_intitial vector and fill with Time in first column and 0s in 2nd

for i = 1:length(Time) %For the length of the Time vector
    Scalar_factor_initial(i) = (neg_d18O(i)-min_d18O)/(max_d18O-min_d18O); %Normalize d18O curve 
end 

%% Part 3: Interpolate Time and scalar factor
%Interpolate time from 143 slices to 661. This will produce ice sheet
%histories every 1000 years. Scalar factor refers to the percentage of ice
%sheet volume that exists for that time slice. A scalar factor of 1 is 100%
%of maximum ice sheet volume (occurs at MIS M2 glacial), whereas a scalar
%factor of 0 indicates that 0% of the maximum ice sheet volume exists at that time
%slice (occurs at MIS MG7 interglacial).

Time_interp = 2950000:1000:3610000; %Interpolate time to every thousand years
x = ((Time*1000).'); %Input Time, put values in row (past to present), and make 100ky scale
y = ((Scalar_factor_initial).'); %Input RSL, put values in row (past to present)
Time_interp = flip(Time_interp); %Flip interpolated time vector from left to right
Scalar_factor_interp = interp1(x,y,Time_interp); %Interpolate scalar factor to full time period

Time = (Time_interp).'; %Flip interpolated time vector from row to column
Scalar_factor = (Scalar_factor_interp).'; %Flip interpolated scalar vector from row to column

Scalar_factor = 1-Scalar_factor; 

%% Part 4: Begin Ice Model Runs for a single GMSL
%% Part 4a: Input ice sheet geometries
%Input MATLAB workspace with individual ice sheet geometries as matrices.
load('Kingetal2024_IceModel_InputIceGeometries1.mat'); %NAIS, EAIS, GrIS, WAIS and EAIS ice geometries
load('Kingetal2024_IceModel_InputIceGeometries1.mat'); %Wilkes Basin, Aurora Basin and Prydz Bay ice geometries
%% Part 4b: Input Ice metrics
%Ice metrics are based off of maximum sea level equivalent (SLE) ice sheet volume at MIS M2 glacial.
%These ice models are run assuming a constant non-marine based ice load in
%Antarctica, therefore BaseAIS_SLE should always be uncommented when
%creating ice histories for EAIS, WAIS, GrIS, EIS, and NAIS. 

BaseAIS_SLE = (Antarctica_IceVolume_DowsettGauss(1)*1e9*(910/1028))/3.635e14;
EAIS_metrics = 14; EAIS_SLE = zeros(661,2);
WAIS_metrics = 5; WAIS_SLE = zeros(661,2);
GrIS_metrics = 7; GrIS_SLE = zeros(661,2);
EIS_metrics = 5; EIS_SLE = zeros(661,2);
NAIS_metrics = 30; NAIS_SLE = zeros(661,2);

%The following are EAIS sub-basin scenarios. These ice models are run assuming a constant non-marine based ice load in
%Antarctica (minus the specidic sub-basin, therefore BaseEAIS*_SLE should always be uncommented when
%creating ice histories for Aurora, Wilkes and Prydz Bay sub-basins. 
BaseEAISAurora_SLE = (EastAntarctica_AuroraBasin_IceVolumeGauss(1)*1e9*(910/1028))/3.635e14;
EAIS_Aurora_metrics = 7.85; EAIS_Aurora_SLE = zeros(661,2);

BaseEAISWilkes_SLE = (EastAntarctica_WilkesBasin_IceVolumeGauss(1)*1e9*(910/1028))/3.635e14;
EAIS_Wilkes_metrics = 5.75; EAIS_Wilkes_SLE = zeros(661,2);

BaseEAISPrydz_SLE = (EastAntarctica_PrydzBay_IceVolumeGauss(1)*1e9*(910/1028))/3.635e14;
EAIS_Prydz_metrics = 2.3; EAIS_Prydz_SLE = zeros(661,2);

% Calcuate SLE proportions of ice for each the 3.61-2.95 Ma time period.
for i = 1:length(Time)
    EAIS_SLE(i,:) = [Time(i),BaseAIS_SLE+(Scalar_factor(i)*EAIS_metrics)];
    EAIS_Aurora_SLE(i,:) = [Time(i),BaseEAISAurora_SLE+(Scalar_factor(i)*EAIS_Aurora_metrics)];
    EAIS_Wilkes_SLE(i,:) = [Time(i),BaseEAISWilkes_SLE+(Scalar_factor(i)*EAIS_Wilkes_metrics)];
    EAIS_Prydz_SLE(i,:) = [Time(i),BaseEAISPrydz_SLE+(Scalar_factor(i)*EAIS_Prydz_metrics)];
    WAIS_SLE(i,:) = [Time(i),(Scalar_factor(i)*WAIS_metrics)];
    GrIS_SLE(i,:) = [Time(i),(Scalar_factor(i)*GrIS_metrics)];
    EIS_SLE(i,:) = [Time(i),(Scalar_factor(i)*EIS_metrics)];
    NAIS_SLE(i,:) = [Time(i),(Scalar_factor(i)*NAIS_metrics)];
end

%% Part 4c: Convert SLE for each ice sheet into an ice volume change in km3
%The SLE changes of each ice are converted to an ice volume in km3 using
%an equation from Goelzer et al., 2020. SLE is the sea equivalent of ice
%above flotation and it's multiplied by the density of the ocean (1028 km
%m^-3) over the density of ice (910 kg m^-3), multiplied by the area of the
%ocean (3.625 x 10^14 m^2) and then divded by 1e9.

EAIS_contribution_km3 = zeros(661,1);
EAIS_Aurora_contribution_km3 = zeros(661,1);
EAIS_Wilkes_contribution_km3 = zeros(661,1);
EAIS_Prydz_contribution_km3 = zeros(661,1);
WAIS_contribution_km3 = zeros(661,1);
GrIS_contribution_km3 = zeros(661,1);
EIS_contribution_km3 = zeros(661,1);
NAIS_contribution_km3 = zeros(661,1);

for i = 1:length(Time)
    EAIS_contribution_km3(i) = (EAIS_SLE(i,2)*(1028/910)*3.635e14)/1e9;
    EAIS_Aurora_contribution_km3(i) = (EAIS_Aurora_SLE(i,2)*(1028/910)*3.635e14)/1e9;
    EAIS_Wilkes_contribution_km3(i) = (EAIS_Wilkes_SLE(i,2)*(1028/910)*3.635e14)/1e9;
    EAIS_Prydz_contribution_km3(i) = (EAIS_Prydz_SLE(i,2)*(1028/910)*3.635e14)/1e9;
    WAIS_contribution_km3(i) = (WAIS_SLE(i,2)*(1028/910)*3.635e14)/1e9;
    GrIS_contribution_km3(i) = (GrIS_SLE(i,2)*(1028/910)*3.635e14)/1e9;
    EIS_contribution_km3(i) = (EIS_SLE(i,2)*(1028/910)*3.635e14)/1e9;
    NAIS_contribution_km3(i) = (NAIS_SLE(i,2)*(1028/910)*3.635e14)/1e9;
end

%% Part 4d: Interpolate Ice Geometries for length of model run
%Ice geometries (geographic extent) for each ice sheet are based off the ~25-40 input Guassian
%Grids, and interpolated for each time slice based on the ice volume in km^3. 

% % Greenland
Greenland_IceGeometries_GaussReshape(isnan(Greenland_IceGeometries_GaussReshape)) = 0; %Repalce NaN values with 0
x = cat(33,Greenland_IceGeometries_GaussReshape(:,:,1),Greenland_IceGeometries_GaussReshape(:,:,2),Greenland_IceGeometries_GaussReshape(:,:,3),Greenland_IceGeometries_GaussReshape(:,:,4),Greenland_IceGeometries_GaussReshape(:,:,5),Greenland_IceGeometries_GaussReshape(:,:,6),Greenland_IceGeometries_GaussReshape(:,:,7),Greenland_IceGeometries_GaussReshape(:,:,8),Greenland_IceGeometries_GaussReshape(:,:,9),Greenland_IceGeometries_GaussReshape(:,:,10),Greenland_IceGeometries_GaussReshape(:,:,11),Greenland_IceGeometries_GaussReshape(:,:,12),Greenland_IceGeometries_GaussReshape(:,:,13),Greenland_IceGeometries_GaussReshape(:,:,14),Greenland_IceGeometries_GaussReshape(:,:,15),Greenland_IceGeometries_GaussReshape(:,:,16),Greenland_IceGeometries_GaussReshape(:,:,17),Greenland_IceGeometries_GaussReshape(:,:,18),Greenland_IceGeometries_GaussReshape(:,:,19),Greenland_IceGeometries_GaussReshape(:,:,20),Greenland_IceGeometries_GaussReshape(:,:,21),Greenland_IceGeometries_GaussReshape(:,:,22),Greenland_IceGeometries_GaussReshape(:,:,23),Greenland_IceGeometries_GaussReshape(:,:,24),Greenland_IceGeometries_GaussReshape(:,:,25),Greenland_IceGeometries_GaussReshape(:,:,26),Greenland_IceGeometries_GaussReshape(:,:,27),Greenland_IceGeometries_GaussReshape(:,:,28),Greenland_IceGeometries_GaussReshape(:,:,29),Greenland_IceGeometries_GaussReshape(:,:,30),Greenland_IceGeometries_GaussReshape(:,:,31),Greenland_IceGeometries_GaussReshape(:,:,32)); %concatenate all ice snapshots
x = permute(x,[33 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32]); %rearrange dimmensions of array
t0 = [Greenland_IceVolumeGauss(:,1) Greenland_IceVolumeGauss(:,2) Greenland_IceVolumeGauss(:,3) Greenland_IceVolumeGauss(:,4) Greenland_IceVolumeGauss(:,5) Greenland_IceVolumeGauss(:,6) Greenland_IceVolumeGauss(:,7) Greenland_IceVolumeGauss(:,8) Greenland_IceVolumeGauss(:,9) Greenland_IceVolumeGauss(:,10) Greenland_IceVolumeGauss(:,11) Greenland_IceVolumeGauss(:,12) Greenland_IceVolumeGauss(:,13) Greenland_IceVolumeGauss(:,14) Greenland_IceVolumeGauss(:,15) Greenland_IceVolumeGauss(:,16) Greenland_IceVolumeGauss(:,17) Greenland_IceVolumeGauss(:,18) Greenland_IceVolumeGauss(:,19) Greenland_IceVolumeGauss(:,20) Greenland_IceVolumeGauss(:,21) Greenland_IceVolumeGauss(:,22) Greenland_IceVolumeGauss(:,23) Greenland_IceVolumeGauss(:,24) Greenland_IceVolumeGauss(:,25) Greenland_IceVolumeGauss(:,26) Greenland_IceVolumeGauss(:,27) Greenland_IceVolumeGauss(:,28) Greenland_IceVolumeGauss(:,29) Greenland_IceVolumeGauss(:,30) Greenland_IceVolumeGauss(:,31) Greenland_IceVolumeGauss(:,32)];

interp_GrIS = []; %Create interp_GrIS matrix
for i = 1:length(Time_interp) %For the length of the Time_interp vector
    interp_GrIS(:,:,i) = interp1(t0,x,GrIS_contribution_km3(i)); %Interpolate each time slice
end

% % North America
NorthAmerica_IceGeometries_GaussReshape(isnan(NorthAmerica_IceGeometries_GaussReshape)) = 0;
x = cat(40,NorthAmerica_IceGeometries_GaussReshape(:,:,1),NorthAmerica_IceGeometries_GaussReshape(:,:,2),NorthAmerica_IceGeometries_GaussReshape(:,:,3),NorthAmerica_IceGeometries_GaussReshape(:,:,4),NorthAmerica_IceGeometries_GaussReshape(:,:,5),NorthAmerica_IceGeometries_GaussReshape(:,:,6),NorthAmerica_IceGeometries_GaussReshape(:,:,7),NorthAmerica_IceGeometries_GaussReshape(:,:,8),NorthAmerica_IceGeometries_GaussReshape(:,:,9),NorthAmerica_IceGeometries_GaussReshape(:,:,10),NorthAmerica_IceGeometries_GaussReshape(:,:,11),NorthAmerica_IceGeometries_GaussReshape(:,:,12),NorthAmerica_IceGeometries_GaussReshape(:,:,13),NorthAmerica_IceGeometries_GaussReshape(:,:,14),NorthAmerica_IceGeometries_GaussReshape(:,:,15),NorthAmerica_IceGeometries_GaussReshape(:,:,16),NorthAmerica_IceGeometries_GaussReshape(:,:,17),NorthAmerica_IceGeometries_GaussReshape(:,:,18),NorthAmerica_IceGeometries_GaussReshape(:,:,19),NorthAmerica_IceGeometries_GaussReshape(:,:,20),NorthAmerica_IceGeometries_GaussReshape(:,:,21),NorthAmerica_IceGeometries_GaussReshape(:,:,22),NorthAmerica_IceGeometries_GaussReshape(:,:,23),NorthAmerica_IceGeometries_GaussReshape(:,:,24),NorthAmerica_IceGeometries_GaussReshape(:,:,25),NorthAmerica_IceGeometries_GaussReshape(:,:,26),NorthAmerica_IceGeometries_GaussReshape(:,:,27),NorthAmerica_IceGeometries_GaussReshape(:,:,28),NorthAmerica_IceGeometries_GaussReshape(:,:,29),NorthAmerica_IceGeometries_GaussReshape(:,:,30),NorthAmerica_IceGeometries_GaussReshape(:,:,31),NorthAmerica_IceGeometries_GaussReshape(:,:,32),NorthAmerica_IceGeometries_GaussReshape(:,:,33),NorthAmerica_IceGeometries_GaussReshape(:,:,34),NorthAmerica_IceGeometries_GaussReshape(:,:,35),NorthAmerica_IceGeometries_GaussReshape(:,:,36),NorthAmerica_IceGeometries_GaussReshape(:,:,37),NorthAmerica_IceGeometries_GaussReshape(:,:,38),NorthAmerica_IceGeometries_GaussReshape(:,:,39));
x = permute(x,[40 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39]);
t0 = [NorthAmerica_IceVolumeGauss(:,1) NorthAmerica_IceVolumeGauss(:,2) NorthAmerica_IceVolumeGauss(:,3) NorthAmerica_IceVolumeGauss(:,4) NorthAmerica_IceVolumeGauss(:,5) NorthAmerica_IceVolumeGauss(:,6) NorthAmerica_IceVolumeGauss(:,7) NorthAmerica_IceVolumeGauss(:,8) NorthAmerica_IceVolumeGauss(:,9) NorthAmerica_IceVolumeGauss(:,10) NorthAmerica_IceVolumeGauss(:,11) NorthAmerica_IceVolumeGauss(:,12) NorthAmerica_IceVolumeGauss(:,13) NorthAmerica_IceVolumeGauss(:,14) NorthAmerica_IceVolumeGauss(:,15) NorthAmerica_IceVolumeGauss(:,16) NorthAmerica_IceVolumeGauss(:,17) NorthAmerica_IceVolumeGauss(:,18) NorthAmerica_IceVolumeGauss(:,19) NorthAmerica_IceVolumeGauss(:,20) NorthAmerica_IceVolumeGauss(:,21) NorthAmerica_IceVolumeGauss(:,22) NorthAmerica_IceVolumeGauss(:,23) NorthAmerica_IceVolumeGauss(:,24) NorthAmerica_IceVolumeGauss(:,25) NorthAmerica_IceVolumeGauss(:,26) NorthAmerica_IceVolumeGauss(:,27) NorthAmerica_IceVolumeGauss(:,28) NorthAmerica_IceVolumeGauss(:,29) NorthAmerica_IceVolumeGauss(:,30) NorthAmerica_IceVolumeGauss(:,31) NorthAmerica_IceVolumeGauss(:,32) NorthAmerica_IceVolumeGauss(:,33) NorthAmerica_IceVolumeGauss(:,34) NorthAmerica_IceVolumeGauss(:,35) NorthAmerica_IceVolumeGauss(:,36) NorthAmerica_IceVolumeGauss(:,37) NorthAmerica_IceVolumeGauss(:,38) NorthAmerica_IceVolumeGauss(:,39)];

interp_NAIS = [];
for i = 1:length(Time_interp)
    interp_NAIS(:,:,i) = interp1(t0,x,NAIS_contribution_km3(i));
end

% %Eurasia
Eurasia_IceGeometries_GaussReshape(isnan(Eurasia_IceGeometries_GaussReshape)) = 0;
x = cat(25,Eurasia_IceGeometries_GaussReshape(:,:,1),Eurasia_IceGeometries_GaussReshape(:,:,2),Eurasia_IceGeometries_GaussReshape(:,:,3),Eurasia_IceGeometries_GaussReshape(:,:,4),Eurasia_IceGeometries_GaussReshape(:,:,5),Eurasia_IceGeometries_GaussReshape(:,:,6),Eurasia_IceGeometries_GaussReshape(:,:,7),Eurasia_IceGeometries_GaussReshape(:,:,8),Eurasia_IceGeometries_GaussReshape(:,:,9),Eurasia_IceGeometries_GaussReshape(:,:,10),Eurasia_IceGeometries_GaussReshape(:,:,11),Eurasia_IceGeometries_GaussReshape(:,:,12),Eurasia_IceGeometries_GaussReshape(:,:,13),Eurasia_IceGeometries_GaussReshape(:,:,14),Eurasia_IceGeometries_GaussReshape(:,:,15),Eurasia_IceGeometries_GaussReshape(:,:,16),Eurasia_IceGeometries_GaussReshape(:,:,17),Eurasia_IceGeometries_GaussReshape(:,:,18),Eurasia_IceGeometries_GaussReshape(:,:,19),Eurasia_IceGeometries_GaussReshape(:,:,20),Eurasia_IceGeometries_GaussReshape(:,:,21),Eurasia_IceGeometries_GaussReshape(:,:,22),Eurasia_IceGeometries_GaussReshape(:,:,23),Eurasia_IceGeometries_GaussReshape(:,:,24));
x = permute(x,[25 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);
t0 = [Eurasia_IceVolumeGauss(:,1) Eurasia_IceVolumeGauss(:,2) Eurasia_IceVolumeGauss(:,3) Eurasia_IceVolumeGauss(:,4) Eurasia_IceVolumeGauss(:,5) Eurasia_IceVolumeGauss(:,6) Eurasia_IceVolumeGauss(:,7) Eurasia_IceVolumeGauss(:,8) Eurasia_IceVolumeGauss(:,9) Eurasia_IceVolumeGauss(:,10) Eurasia_IceVolumeGauss(:,11) Eurasia_IceVolumeGauss(:,12) Eurasia_IceVolumeGauss(:,13) Eurasia_IceVolumeGauss(:,14) Eurasia_IceVolumeGauss(:,15) Eurasia_IceVolumeGauss(:,16) Eurasia_IceVolumeGauss(:,17) Eurasia_IceVolumeGauss(:,18) Eurasia_IceVolumeGauss(:,19) Eurasia_IceVolumeGauss(:,20) Eurasia_IceVolumeGauss(:,21) Eurasia_IceVolumeGauss(:,22) Eurasia_IceVolumeGauss(:,23) Eurasia_IceVolumeGauss(:,24)];

interp_EIS = [];
for i = 1:length(Time_interp)
    interp_EIS(:,:,i) = interp1(t0,x,EIS_contribution_km3(i));
end

% %East Antarctica
EastAntarctica_IceGeometries_GaussReshape(isnan(EastAntarctica_IceGeometries_GaussReshape)) = 0;
x = cat(26,EastAntarctica_IceGeometries_GaussReshape(:,:,1),EastAntarctica_IceGeometries_GaussReshape(:,:,2),EastAntarctica_IceGeometries_GaussReshape(:,:,3),EastAntarctica_IceGeometries_GaussReshape(:,:,4),EastAntarctica_IceGeometries_GaussReshape(:,:,5),EastAntarctica_IceGeometries_GaussReshape(:,:,6),EastAntarctica_IceGeometries_GaussReshape(:,:,7),EastAntarctica_IceGeometries_GaussReshape(:,:,8),EastAntarctica_IceGeometries_GaussReshape(:,:,9),EastAntarctica_IceGeometries_GaussReshape(:,:,10),EastAntarctica_IceGeometries_GaussReshape(:,:,11),EastAntarctica_IceGeometries_GaussReshape(:,:,12),EastAntarctica_IceGeometries_GaussReshape(:,:,13),EastAntarctica_IceGeometries_GaussReshape(:,:,14),EastAntarctica_IceGeometries_GaussReshape(:,:,15),EastAntarctica_IceGeometries_GaussReshape(:,:,16),EastAntarctica_IceGeometries_GaussReshape(:,:,17),EastAntarctica_IceGeometries_GaussReshape(:,:,18),EastAntarctica_IceGeometries_GaussReshape(:,:,19),EastAntarctica_IceGeometries_GaussReshape(:,:,20),EastAntarctica_IceGeometries_GaussReshape(:,:,21),EastAntarctica_IceGeometries_GaussReshape(:,:,22),EastAntarctica_IceGeometries_GaussReshape(:,:,23),EastAntarctica_IceGeometries_GaussReshape(:,:,24),EastAntarctica_IceGeometries_GaussReshape(:,:,25));
x = permute(x,[26 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25]);
t0 = [EastAntarctica_IceVolumeGauss(:,1) EastAntarctica_IceVolumeGauss(:,2) EastAntarctica_IceVolumeGauss(:,3) EastAntarctica_IceVolumeGauss(:,4) EastAntarctica_IceVolumeGauss(:,5) EastAntarctica_IceVolumeGauss(:,6) EastAntarctica_IceVolumeGauss(:,7) EastAntarctica_IceVolumeGauss(:,8) EastAntarctica_IceVolumeGauss(:,9) EastAntarctica_IceVolumeGauss(:,10) EastAntarctica_IceVolumeGauss(:,11) EastAntarctica_IceVolumeGauss(:,12) EastAntarctica_IceVolumeGauss(:,13) EastAntarctica_IceVolumeGauss(:,14) EastAntarctica_IceVolumeGauss(:,15) EastAntarctica_IceVolumeGauss(:,16) EastAntarctica_IceVolumeGauss(:,17) EastAntarctica_IceVolumeGauss(:,18) EastAntarctica_IceVolumeGauss(:,19) EastAntarctica_IceVolumeGauss(:,20) EastAntarctica_IceVolumeGauss(:,21) EastAntarctica_IceVolumeGauss(:,22) EastAntarctica_IceVolumeGauss(:,23) EastAntarctica_IceVolumeGauss(:,24) EastAntarctica_IceVolumeGauss(:,25)];

interp_EAIS = [];
for i = 1:length(Time_interp)
    interp_EAIS(:,:,i) = interp1(t0,x,EAIS_contribution_km3(i));
end

% %Aurora Basin
EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(isnan(EastAntarctica_AuroraBasin_IceGeometries_GaussReshape)) = 0;
x = cat(10,EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,1),EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,2),EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,3),EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,4),EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,5),EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,6),EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,7),EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,8),EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,9));
x = permute(x,[10 1 2 3 4 5 6 7 8 9]);
t0 = [EastAntarctica_AuroraBasin_IceVolumeGauss(:,1) EastAntarctica_AuroraBasin_IceVolumeGauss(:,2) EastAntarctica_AuroraBasin_IceVolumeGauss(:,3) EastAntarctica_AuroraBasin_IceVolumeGauss(:,4) EastAntarctica_AuroraBasin_IceVolumeGauss(:,5) EastAntarctica_AuroraBasin_IceVolumeGauss(:,6) EastAntarctica_AuroraBasin_IceVolumeGauss(:,7) EastAntarctica_AuroraBasin_IceVolumeGauss(:,8) EastAntarctica_AuroraBasin_IceVolumeGauss(:,9)];

interp_EAIS_Aurora = [];
for i = 1:length(Time_interp)
    interp_EAIS_Aurora(:,:,i) = interp1(t0,x,EAIS_Aurora_contribution_km3(i));
end

% %Wilkes Basin
EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(isnan(EastAntarctica_WilkesBasin_IceGeometries_GaussReshape)) = 0;
x = cat(9,EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(:,:,1),EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(:,:,2),EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(:,:,3),EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(:,:,4),EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(:,:,5),EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(:,:,6),EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(:,:,7),EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(:,:,8));
x = permute(x,[9 1 2 3 4 5 6 7 8]);
t0 = [EastAntarctica_WilkesBasin_IceVolumeGauss(:,1) EastAntarctica_WilkesBasin_IceVolumeGauss(:,2) EastAntarctica_WilkesBasin_IceVolumeGauss(:,3) EastAntarctica_WilkesBasin_IceVolumeGauss(:,4) EastAntarctica_WilkesBasin_IceVolumeGauss(:,5) EastAntarctica_WilkesBasin_IceVolumeGauss(:,6) EastAntarctica_WilkesBasin_IceVolumeGauss(:,7) EastAntarctica_WilkesBasin_IceVolumeGauss(:,8)];

interp_EAIS_Wilkes = [];
for i = 1:length(Time_interp)
    interp_EAIS_Wilkes(:,:,i) = interp1(t0,x,EAIS_Wilkes_contribution_km3(i));
end

% %Prydz Basin
EastAntarctica_PrydzBay_IceGeometries_GaussReshape(isnan(EastAntarctica_PrydzBay_IceGeometries_GaussReshape)) = 0;
x = cat(7,EastAntarctica_PrydzBay_IceGeometries_GaussReshape(:,:,1),EastAntarctica_PrydzBay_IceGeometries_GaussReshape(:,:,2),EastAntarctica_PrydzBay_IceGeometries_GaussReshape(:,:,3),EastAntarctica_PrydzBay_IceGeometries_GaussReshape(:,:,4),EastAntarctica_PrydzBay_IceGeometries_GaussReshape(:,:,5),EastAntarctica_PrydzBay_IceGeometries_GaussReshape(:,:,6));
x = permute(x,[7 1 2 3 4 5 6]);
t0 = [EastAntarctica_PrydzBay_IceVolumeGauss(:,1) EastAntarctica_PrydzBay_IceVolumeGauss(:,2) EastAntarctica_PrydzBay_IceVolumeGauss(:,3) EastAntarctica_PrydzBay_IceVolumeGauss(:,4) EastAntarctica_PrydzBay_IceVolumeGauss(:,5) EastAntarctica_PrydzBay_IceVolumeGauss(:,6)];

interp_EAIS_Prydz = [];
for i = 1:length(Time_interp)
    interp_EAIS_Prydz(:,:,i) = interp1(t0,x,EAIS_Prydz_contribution_km3(i));
end

% %West Antarctica
WestAntarctica_IceGeometries_GaussReshape(isnan(WestAntarctica_IceGeometries_GaussReshape)) = 0;
x = cat(26,WestAntarctica_IceGeometries_GaussReshape(:,:,1),WestAntarctica_IceGeometries_GaussReshape(:,:,2),WestAntarctica_IceGeometries_GaussReshape(:,:,3),WestAntarctica_IceGeometries_GaussReshape(:,:,4),WestAntarctica_IceGeometries_GaussReshape(:,:,5),WestAntarctica_IceGeometries_GaussReshape(:,:,6),WestAntarctica_IceGeometries_GaussReshape(:,:,7),WestAntarctica_IceGeometries_GaussReshape(:,:,8),WestAntarctica_IceGeometries_GaussReshape(:,:,9),WestAntarctica_IceGeometries_GaussReshape(:,:,10),WestAntarctica_IceGeometries_GaussReshape(:,:,11),WestAntarctica_IceGeometries_GaussReshape(:,:,12),WestAntarctica_IceGeometries_GaussReshape(:,:,13),WestAntarctica_IceGeometries_GaussReshape(:,:,14),WestAntarctica_IceGeometries_GaussReshape(:,:,15),WestAntarctica_IceGeometries_GaussReshape(:,:,16),WestAntarctica_IceGeometries_GaussReshape(:,:,17),WestAntarctica_IceGeometries_GaussReshape(:,:,18),WestAntarctica_IceGeometries_GaussReshape(:,:,19),WestAntarctica_IceGeometries_GaussReshape(:,:,20),WestAntarctica_IceGeometries_GaussReshape(:,:,21),WestAntarctica_IceGeometries_GaussReshape(:,:,22),WestAntarctica_IceGeometries_GaussReshape(:,:,23),WestAntarctica_IceGeometries_GaussReshape(:,:,24),WestAntarctica_IceGeometries_GaussReshape(:,:,25));
x = permute(x,[26 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25]);
t0 = [WestAntarctica_IceVolumeGauss(:,1) WestAntarctica_IceVolumeGauss(:,2) WestAntarctica_IceVolumeGauss(:,3) WestAntarctica_IceVolumeGauss(:,4) WestAntarctica_IceVolumeGauss(:,5) WestAntarctica_IceVolumeGauss(:,6) WestAntarctica_IceVolumeGauss(:,7) WestAntarctica_IceVolumeGauss(:,8) WestAntarctica_IceVolumeGauss(:,9) WestAntarctica_IceVolumeGauss(:,10) WestAntarctica_IceVolumeGauss(:,11) WestAntarctica_IceVolumeGauss(:,12) WestAntarctica_IceVolumeGauss(:,13) WestAntarctica_IceVolumeGauss(:,14) WestAntarctica_IceVolumeGauss(:,15) WestAntarctica_IceVolumeGauss(:,16) WestAntarctica_IceVolumeGauss(:,17) WestAntarctica_IceVolumeGauss(:,18) WestAntarctica_IceVolumeGauss(:,19) WestAntarctica_IceVolumeGauss(:,20) WestAntarctica_IceVolumeGauss(:,21) WestAntarctica_IceVolumeGauss(:,22) WestAntarctica_IceVolumeGauss(:,23) WestAntarctica_IceVolumeGauss(:,24) WestAntarctica_IceVolumeGauss(:,25)];

interp_WAIS = [];
for i = 1:length(Time_interp)
    interp_WAIS(:,:,i) = interp1(t0,x,WAIS_contribution_km3(i));
end

%% Part 5: Combine everything onto Global Grid

%Create base Antarctica matrices
BaseAntarctica_IceGeometries_GaussReshape = zeros(512,256,661); %Create BaseAntarctica_IceGeometries_GaussReshape file
for i = 1:661
    BaseAntarctica_IceGeometries_GaussReshape(:,:,i) = Antarctica_IceGeometries_DowsettGaussReshape(:,:,1); %Repalce every matrix with the original base Antarctic ice geometry
end
interp_EAIS_Wilkes(:,:,41)=EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(:,:,1); %At peak interglacial periods (where there is no ice left in sub-basin) the ice sheet geometry should match the base Antarctic ice geometry
interp_EAIS_Wilkes(:,:,131)=EastAntarctica_WilkesBasin_IceGeometries_GaussReshape(:,:,1);

interp_EAIS_Aurora(:,:,41)=EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,1);
interp_EAIS_Aurora(:,:,131)=EastAntarctica_AuroraBasin_IceGeometries_GaussReshape(:,:,1);

interp_EAIS_Prydz(:,:,41)=EastAntarctica_PrydzBay_IceGeometries_GaussReshape(:,:,1);
interp_EAIS_Prydz(:,:,131)=EastAntarctica_PrydzBay_IceGeometries_GaussReshape(:,:,1);

%Create matrices with the base Antarctica geometry and interpolated individual ice
%sheet geometries over the 3.61-2.95 time period
GaussGrid_AllIce_NAIS =  BaseAntarctica_IceGeometries_GaussReshape + interp_NAIS;
GaussGrid_AllIce_GrIS =  BaseAntarctica_IceGeometries_GaussReshape + interp_GrIS; 
GaussGrid_AllIce_EAIS =  BaseAntarctica_IceGeometries_GaussReshape + interp_EAIS; 
GaussGrid_AllIce_WAIS =  BaseAntarctica_IceGeometries_GaussReshape + interp_WAIS; 
GaussGrid_AllIce_EIS =  BaseAntarctica_IceGeometries_GaussReshape + interp_EIS; 
GaussGrid_AllIce_EAIS_Wilkes = interp_EAIS_Wilkes;
GaussGrid_AllIce_EAIS_Aurora = interp_EAIS_Aurora;
GaussGrid_AllIce_EAIS_Prydz = interp_EAIS_Prydz;
%% Part 6: Create binary ice files from Gauss Matrix
%This part of the code is set up to run assuming that there is already a
%"times" binary file and a series of "ice" binary files within an ice sheet folder. I
%would suggest making a folder for each ice sheet and prepopulating it with
%those files, then the code will overwrite them with the new data.

%The sea level code needs the time files structured in years before present
Time_interp = 2950000:1000:3610000; %Create time vector from 2.95 - 3.61 Ma, with a time slive every 1000 years
Time_interp = flip(Time_interp); %Flip vector so that it starts with 3.61 Ma
Time = -(Time_interp).'; %Convert Time_interp to column instead of row and add a negative to each time slice

%NAIS
fileID = fopen('.../NAIS/times','w'); %Open "times" file
fwrite(fileID,Time,'float32'); %Overwrite file with data from Times
fclose(fileID); %Close file

for i = 1:661 %For each of the 661 time slices
    formatSpec = '.../NAIS/ice%d'; %Format "ice" file with number of time slice (e.g., for first loop file should be "ice1")
    filename = sprintf(formatSpec,i);
    fileID = fopen(filename,'w'); %Open "ice" file
    fwrite(fileID,GaussGrid_AllIce_NAIS(:,:,i),'float32'); %Overwrite file with data from GaussGrid_AllIce_NAIS
    fclose(fileID); %Close file
end

%GrIS
fileID = fopen('/Volumes/Extreme SSD/IceModelOutput/GrIS/times','w');
fwrite(fileID,Time,'float32');
fclose(fileID);

for i = 1:661
    formatSpec = '/Volumes/Extreme SSD/IceModelOutput/GrIS/ice%d';
    filename = sprintf(formatSpec,i);
    fileID = fopen(filename,'w');
    fwrite(fileID,GaussGrid_AllIce_GrIS(:,:,i),'float32');
    fclose(fileID);
end

%EAIS
fileID = fopen('/Volumes/Extreme SSD/IceModelOutput/EAIS/times','w');
fwrite(fileID,Time,'float32');
fclose(fileID);

for i = 1:661
    formatSpec = '/Volumes/Extreme SSD/IceModelOutput/EAIS/ice%d';
    filename = sprintf(formatSpec,i);
    fileID = fopen(filename,'w');
    fwrite(fileID,GaussGrid_AllIce_EAIS(:,:,i),'float32');
    fclose(fileID);
end

%WAIS
fileID = fopen('/Volumes/Extreme SSD/IceModelOutput/WAIS/times','w');
fwrite(fileID,Time,'float32');
fclose(fileID);

for i = 1:661
    formatSpec = '/Volumes/Extreme SSD/IceModelOutput/WAIS/ice%d';
    filename = sprintf(formatSpec,i);
    fileID = fopen(filename,'w');
    fwrite(fileID,GaussGrid_AllIce_WAIS(:,:,i),'float32');
    fclose(fileID);
end

%EIS
fileID = fopen('/Volumes/Extreme SSD/IceModelOutput/EIS/times','w');
fwrite(fileID,Time,'float32');
fclose(fileID);

for i = 1:661
    formatSpec = '/Volumes/Extreme SSD/IceModelOutput/EIS/ice%d';
    filename = sprintf(formatSpec,i);
    fileID = fopen(filename,'w');
    fwrite(fileID,GaussGrid_AllIce_EIS(:,:,i),'float32');
    fclose(fileID);
end

%Wilkes Basin
fileID = fopen('/Volumes/Extreme SSD/IceModelOutput/WilkesBasin/times','w');
fwrite(fileID,Time,'float32');
fclose(fileID);

for i = 1:661
    formatSpec = '/Volumes/Extreme SSD/IceModelOutput/WilkesBasin/ice%d';
    filename = sprintf(formatSpec,i);
    fileID = fopen(filename,'w');
    fwrite(fileID,GaussGrid_AllIce_EAIS_Wilkes(:,:,i),'float32');
    fclose(fileID);
end

%Aurora Basin
fileID = fopen('/Volumes/Extreme SSD/IceModelOutput/AuroraBasin/times','w');
fwrite(fileID,Time,'float32');
fclose(fileID);

for i = 1:661
    formatSpec = '/Volumes/Extreme SSD/IceModelOutput/AuroraBasin/ice%d';
    filename = sprintf(formatSpec,i);
    fileID = fopen(filename,'w');
    fwrite(fileID,GaussGrid_AllIce_EAIS_Aurora(:,:,i),'float32');
    fclose(fileID);
end

%Prydz Bay
fileID = fopen('/Volumes/Extreme SSD/IceModelOutput/PrydzBay/times','w');
fwrite(fileID,Time,'float32');
fclose(fileID);

for i = 1:661
    formatSpec = '/Volumes/Extreme SSD/IceModelOutput/PrydzBay/ice%d';
    filename = sprintf(formatSpec,i);
    fileID = fopen(filename,'w');
    fwrite(fileID,GaussGrid_AllIce_EAIS_Prydz(:,:,i),'float32');
    fclose(fileID);
end