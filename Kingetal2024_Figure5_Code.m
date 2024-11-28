%% Figure 5
%Code to plot Figure 5 in King et al., 2024. Inputs are the "Summary" Sea
%Level Histories from all 27 earth models outputted by the Sea Level Model.
%The Summary files are text files with the GMSLP values from each earth
%model run (in the order that they appear in Table 1).
%The code outputs 3 separate figures of box and whisker plots that were
%combined in Adobe Illustrator to create Figure 5 in the manuscript.

%% Input SUMMARY folder of sea level histories from all 27 earth model runs
mainfolder =  '.../SUMMARY_Output'; %Input SUMMARY folder
filelist_SUMMARY_Output = dir(fullfile(mainfolder, '**')); %Create a list of all the files within the SUMMARY folder

for i = 4:220 %length of items the SUMMARY folder
    rootfolder = getfield(filelist_SUMMARY_Output(i),'folder'); %Extract the name of the file as a string
    filename = fullfile(rootfolder,filelist_SUMMARY_Output(i).name); %Create filename vector with string
    dataMatrix = readmatrix(filename); %Import the matrix data from the file
    
    %Rename the file with the name of the model run
    out = regexp(filename,'/','split');
    outstring = string(out(6));
    match = ["%","."];
    newStr = erase(outstring,match);
    formatSpec = '%s_12';
    varname5 = sprintf(formatSpec,newStr);
    
    if length(dataMatrix) == 9
        dataMatrix(1)=[];
    else
    end
    
    %Create separate vectors for each GMSLP, Enewetak Atoll, New Zealand,
    %and Virginia value
    eval(sprintf('%s_GMSLP = %d',newStr,dataMatrix(:,1)));
    eval(sprintf('%s_Enewetak = %d',newStr,dataMatrix(:,4)));
    eval(sprintf('%s_NewZealand = %d',newStr,dataMatrix(:,6)));
    eval(sprintf('%s_Virginia = %d',newStr,dataMatrix(:,8)));
end
%% Extract GMSLP values from text files

%Create vectors with the order of values for the earth models. 
EarthModels3 = [72;72;72;72;72;72;72;72;72;96;96;96;96;96;96;96;96;96;125;125;125;125;125;125;125;125;125]; %Lithospheric thickness values
EarthModels1 = [3;3;3;5;5;5;8;8;8;3;3;3;5;5;5;8;8;8;3;3;3;5;5;5;8;8;8]; %Upper mantle viscosity values
EarthModels2 = [5;10;30;5;10;30;5;10;30;5;10;30;5;10;30;5;10;30;5;10;30;5;10;30;5;10;30]; %Lower mantle viscosity values

%% Create Normalized Vectors by each geographic location

%New Zealand
NormalizedNorthAmerica_NewZealand = []; %Create empty vector containing GMSLP values all 27 earth model runs for each ice sheet 
NormalizedEastAntarctica_NewZealand = [];
NormalizedGreenland_NewZealand = [];
NormalizedWestAntarctica_NewZealand = [];
NormalizedEurasia_NewZealand = [];
NormalizedAuroraBasin_NewZealand = [];
NormalizedWilkesBasin_NewZealand = [];
NormalizedPrydzBay_NewZealand = [];

%Populate those vectors with the values from the corresponding ice sheet
%and geographic location
for i = 1:length(EarthModels1)
    eval(sprintf('NormalizedNorthAmerica_NewZealand=[NormalizedNorthAmerica_NewZealand;output_summary_%dCp%d%d_NorthAmerica35mtxt_NewZealand]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedEastAntarctica_NewZealand=[NormalizedEastAntarctica_NewZealand;output_summary_%dCp%d%d_EastAntarcticatxt_NewZealand]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedGreenland_NewZealand=[NormalizedGreenland_NewZealand;output_summary_%dCp%d%d_Greenlandtxt_NewZealand]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedWestAntarctica_NewZealand=[NormalizedWestAntarctica_NewZealand;output_summary_%dCp%d%d_WestAntarcticatxt_NewZealand]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedEurasia_NewZealand=[NormalizedEurasia_NewZealand;output_summary_%dCp%d%d_Eurasiatxt_NewZealand]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedAuroraBasin_NewZealand=[NormalizedAuroraBasin_NewZealand;output_summary_%dCp%d%d_AuroraBasintxt_NewZealand]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedWilkesBasin_NewZealand=[NormalizedWilkesBasin_NewZealand;output_summary_%dCp%d%d_WilkesBasintxt_NewZealand]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedPrydzBay_NewZealand=[NormalizedPrydzBay_NewZealand;output_summary_%dCp%d%d_PrydzBaytxt_NewZealand]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));

end
   

% Virginia 
NormalizedNorthAmerica_Virginia = [];
NormalizedEastAntarctica_Virginia = [];
NormalizedGreenland_Virginia = [];
NormalizedWestAntarctica_Virginia = [];
NormalizedEurasia_Virginia = [];
NormalizedAuroraBasin_Virginia = [];
NormalizedWilkesBasin_Virginia = [];
NormalizedPrydzBay_Virginia = [];

for i = 1:length(EarthModels1)
    eval(sprintf('NormalizedNorthAmerica_Virginia=[NormalizedNorthAmerica_Virginia;output_summary_%dCp%d%d_NorthAmerica35mtxt_Virginia]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedEastAntarctica_Virginia=[NormalizedEastAntarctica_Virginia;output_summary_%dCp%d%d_EastAntarcticatxt_Virginia]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedGreenland_Virginia=[NormalizedGreenland_Virginia;output_summary_%dCp%d%d_Greenlandtxt_Virginia]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedWestAntarctica_Virginia=[NormalizedWestAntarctica_Virginia;output_summary_%dCp%d%d_WestAntarcticatxt_Virginia]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedEurasia_Virginia=[NormalizedEurasia_Virginia;output_summary_%dCp%d%d_Eurasiatxt_Virginia]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedAuroraBasin_Virginia=[NormalizedAuroraBasin_Virginia;output_summary_%dCp%d%d_AuroraBasintxt_Virginia]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedWilkesBasin_Virginia=[NormalizedWilkesBasin_Virginia;output_summary_%dCp%d%d_WilkesBasintxt_Virginia]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedPrydzBay_Virginia=[NormalizedPrydzBay_Virginia;output_summary_%dCp%d%d_PrydzBaytxt_Virginia]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));

end


% Enewetak Atoll
NormalizedNorthAmerica_Enewetak = [];
NormalizedEastAntarctica_Enewetak = [];
NormalizedGreenland_Enewetak = [];
NormalizedWestAntarctica_Enewetak = [];
NormalizedEurasia_Enewetak = [];
NormalizedAuroraBasin_Enewetak = [];
NormalizedWilkesBasin_Enewetak = [];
NormalizedPrydzBay_Enewetak = [];

for i = 1:length(EarthModels1)
    eval(sprintf('NormalizedNorthAmerica_Enewetak=[NormalizedNorthAmerica_Enewetak;output_summary_%dCp%d%d_NorthAmerica35mtxt_Enewetak]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedEastAntarctica_Enewetak=[NormalizedEastAntarctica_Enewetak;output_summary_%dCp%d%d_EastAntarcticatxt_Enewetak]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedGreenland_Enewetak=[NormalizedGreenland_Enewetak;output_summary_%dCp%d%d_Greenlandtxt_Enewetak]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedWestAntarctica_Enewetak=[NormalizedWestAntarctica_Enewetak;output_summary_%dCp%d%d_WestAntarcticatxt_Enewetak]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedEurasia_Enewetak=[NormalizedEurasia_Enewetak;output_summary_%dCp%d%d_Eurasiatxt_Enewetak]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedAuroraBasin_Enewetak=[NormalizedAuroraBasin_Enewetak;output_summary_%dCp%d%d_AuroraBasintxt_Enewetak]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedWilkesBasin_Enewetak=[NormalizedWilkesBasin_Enewetak;output_summary_%dCp%d%d_WilkesBasintxt_Enewetak]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));
    eval(sprintf('NormalizedPrydzBay_Enewetak=[NormalizedPrydzBay_Enewetak;output_summary_%dCp%d%d_PrydzBaytxt_Enewetak]',EarthModels3(i),EarthModels1(i),EarthModels2(i)));

end


%% Create box and whisker plots 
%Box and whisker plots are created for the variation in GMSLP for each ice 
%sheet for each geographic location

%Combine all the data for each geographic location
VirginiaPlot = [NormalizedNorthAmerica_Virginia,NormalizedGreenland_Virginia,NormalizedEurasia_Virginia,NormalizedEastAntarctica_Virginia,NormalizedWestAntarctica_Virginia,NormalizedAuroraBasin_Virginia,NormalizedWilkesBasin_Virginia,NormalizedPrydzBay_Virginia];
EnewetakPlot = [NormalizedNorthAmerica_Enewetak,NormalizedGreenland_Enewetak,NormalizedEurasia_Enewetak,NormalizedEastAntarctica_Enewetak,NormalizedWestAntarctica_Enewetak,NormalizedAuroraBasin_Enewetak,NormalizedWilkesBasin_Enewetak,NormalizedPrydzBay_Enewetak];
NewZealandPlot = [NormalizedNorthAmerica_NewZealand,NormalizedGreenland_NewZealand,NormalizedEurasia_NewZealand,NormalizedEastAntarctica_NewZealand,NormalizedWestAntarctica_NewZealand,NormalizedAuroraBasin_NewZealand,NormalizedWilkesBasin_NewZealand,NormalizedPrydzBay_NewZealand];

%Plot all ice sheets and earth model runs for New Zealand
figure
boxplot(NewZealandPlot)
ylim([.25 1.2])
title("New Zealand",'FontSize',18)
xticklabels({'North America','Greenland','Eurasia','East Antarcica','West Antarctica','Aurora Basin','Wilkes Basin','Prydz Bay'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)

%Plot all ice sheets and earth model runs for Virginia
figure
boxplot(VirginiaPlot)
ylim([.25 1.2])
title("Virginia",'FontSize',18)
xticklabels({'North America','Greenland','Eurasia','East Antarcica','West Antarctica','Aurora Basin','Wilkes Basin','Prydz Bay'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)

%Plot all ice sheets and earth model runs for Enewetak Atoll
figure
boxplot(EnewetakPlot)
ylim([.25 1.2])
title("Enewetak Atoll",'FontSize',18)
xticklabels({'North America','Greenland','Eurasia','East Antarcica','West Antarctica','Aurora Basin','Wilkes Basin','Prydz Bay'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)