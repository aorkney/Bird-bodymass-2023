# Bird-bodymass-2023
This repository contains a collection of scripts that investigate how body mass variety in birds structures the modular evolutionary organisation of the avian skeleton, and the strength of its ecological adaptation to flight style variety.

The scripts herein, which were all run in R verson 4.2.1, reproduce the analysis in the submitted manuscript Orkney & Hedrick 2023; 'Penetrating the Palimpsest; bodymass drives divergent patterns of modular evolution in the avian wing and trunk,'

The script 'Compile_datasets_22_10_2022_clean.R' will draw upon supplied or cited source files to generate and save datasets for analysis in subsequent scripts.

Figure 1:
Figure 1 is a collection of line plots which will help the user visualise how the strength of evolutionary integration between bones within different regions/modules of the avian skeleton change as avian body mass increases. This is achieved by generating many taxonomic subsets of birds with different mean body mass values, and repeating analyses to determin evolutionary integration.

a) Run script 'Residual_Csize_10_08_2022.R' to generate a dataset of avian bone sizes (centroid size computed from landmark constellations), corrected for patterns of allometric scaling (size-dependent scaling; this must be removed because the role of body mass is being investigated and size-dependent scaling could obfuscate such patterns).

b) Run script 'Integration_v_body_mass_11_08_2022.R' to run the analysis.

c) Run script 'Plot_integration_v_body_mass_14_11_2022.R' to plot the analysis output (Figure 1).


Figure 2:
Figure 2 is a collection of raster images which demonstrate whether evolutionary integration within and between the wing and trunk skeletal modules changes significantly as a function of body mass. This is achieved by producing many taxonomic subsets of birds with different mass ranges and asking whether evolutionary covariance between selected combinations of bones is significantly stronger or weaker in lighter or more massive taxonomic subsets. 

a) Run script 'pairwise_sig_compare_28_12_2022.R' to undertake analysis and plot Figure 2.

Figure 3: 
Figure 3 employs an alternative mathematical approach to test whether evolutionary integration changes as a function of body mass. In this analysis the evolutionary covariances between bones within the wing (carpometacarpus and humerus), between the size of the carpometacarpus and flight style ecology, and between the 3-D shape of the carpometacarpus and flight style ecology, are computed. The major axes of covariance are extracted and the Euclidean distance of individual species to these axes are found. Breusch-pagan tests for heteroskedasticity are then conducted to determine whether the scatter from the major axis of covariance depends on body mass (stronger/weaker integration with higher body mass will result in a narrower/wider scatter).
The results are represented as a series of decorated scatter plots. 

a) Run script 'Wing_plot_04_12_2022.R' to undertake analysis and plot Figure 3.


Figure 4:
Figure 4 employs the same approach as Figure 3, but in this instance evolutionary covariances between bones within the trunk (sternum and scapula), between the trunk and wing (sternum and carpometacarpus) and between 3-D shape of the coracoid (a trunk bone the shape of which covaries with flight style) and flight style ecology. 

a) Run script 'Trunk_plot_04_11_2022.R' to undertake analysis and plot Figure 4.
