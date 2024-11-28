# plioceneSeaLevel
Data and code associated with "The Geometry of Sea-Level Change Across a Mid Pliocene Glacial Cycle"

# Figures
Contains MATLAB .m files for Figures 2-7. 
Figure 1 was entirely drawn in Adobe Illustrator so there is no code associated with it. Figures 3,4,6 and 7 will need the cbrewer MATLAB package in order to plot the Guassian Grids.

# IceModel
Contains MATLAB .m, .mat files as well as a .txt and several binary files. 
_Kingetal2024_IceModel_Code.m_ is the code to run the ice model. 
_Kingetal2024_IceModel_Inputs18O.txt_ is an input into the ice model. It is a .txt file containing a portion of the LR04 d18O stack. 
_Kingetal2024_IceModel_InputIceGeometries1.m_ and _Kingetal2024_IceModel_InputIceGeometries1.m_ are inputs into the ice model. They are MATLAB workspaces containing the ice snapshot geometries.
_Kingetal2024_IceModel_OutputIceHistories.zip_ is the output from the ice model. It is a zipped folder containing the ice316 and ice456 binary files from each ice sheet. The full ice histories are too large to include here, but are avaialable upon request.

# SeaLevelModel
Contains a FORTRAN .f90 file as well as several binary and .txt files.
_Kingetal2024_SeaLevelModel_Code.f90_ is the code to run the sea level model. It was written by Sam Goldberg (University of Miami).
_Kingetal2024_SeaLevelModel_InputEarthModels.zip_ is an input into the sea level model. It is a zipped folder containing the earth model binary files. 
_Kingetal2024_SeaLevelModel_OutputSLGrids.zip_ is an output from the sea level model. It is a zipped folder containing binary sea level grids for each ice sheet scenario. The grids are the the difference (i.e., sea-level change) between the KM3 interglacial and M2 glacial. One set of grids are the relative sea-level change, while the other are the normalized sea-level change.
_Kingetal2024_SeaLevelModel_OutputSummaryFiles.zip_ is an output from the sea level model. It is a zipped folder containing text files with the GMSLP values from each earth model run (in the order that they appear in Table 1).
*Note that the other input into the sea level model are the output ice histories produced by the ice model.

