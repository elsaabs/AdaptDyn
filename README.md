This repository contains the data needed to run the code files, and the code files to analyze the data and get the plots in the manuscript. 

Data files
+ EFC_Carbohydrate_Lignin_Temp.csv: omics data table. 
+ enzyme_activity_combined.xlsx: Enzyme assay data table.
+ /nc_files/cLeaf_Lmon_CCSM4_rcp85_r1i1p1_200601-210012.nc: CCSM4 RCP8.5 temperature and litter predictions. 

Code files
+ Rfile 27: Gives the predicted global map of soil temperature from 2010 to 2100 (fig. S1, fig. S4). 
+ Rfile 32: Gives the global average soil temperature as a function of time between 2010 to 2100 (Fig. 2, temperature curve). 
+ Rfile 33: Gives the predicted change in global litter per decade between 2010 and 2100 (fig. S7). 
+ Rfile 44: Gives the global maps of 14 to 5 biomes (fig. S3).
+ Rfile 108: Gives the global map of adapted enzyme allocation (phiecoevo) in 2010, from which I derive the global map of stable points in 2010 (Fig. 4a).
+ Rfile 109: Checks that the global map of stable points in 2010 with eco-evolution (phi varies with temperature) and without evolution (phi = phiecoevo) are the same (as we start with the same initial conditions in 2010). 
+ Rfile 110: Gives the global maps of optimal enzyme allocation and of SOC with eco-evolution each decade between 2010 and 2100 (Fig. 4c, Fig. 3b, fig. S2).
+ Rfile 111: Gives the global map of SOC without eco-evolution each decade between 2010 and 2100 (Fig. 3a).
+ Rfile 112: Gives the total global SOC as a function of time in years between 2010 and 2100 with vs without eco-evolution (Fig. 2). 
+ Code_for_Science_Advances.nb: Mathematica code file for model analysis (Fig. 3c, Fig. 4b,d, fig. S6). 
+ Weintraub_data_analysis_for_ms.ipynb: jupyter notebook for enzyme assay data analysis (fig. S11a-c, fig. S13). 
+ Yang_data_analysis.ipynb: jupyter notebook for omics data analysis (fig. S10, fig. S11d-f). 
