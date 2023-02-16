# Bird-bodymass-2023
This repository contains a collection of scripts that investigate how body mass variety in birds structures the modular evolutionary organisation of the avian skeleton, and the strength of its ecological adaptation to flight style variety.

The scripts herein, which were all run in R verson 4.2.1, reproduce the analysis in the submitted manuscript Orkney & Hedrick 2023; 'Penetrating the Palimpsest; bodymass drives divergent patterns of modular evolution in the avian wing and trunk,'

The script 'Compile_datasets_22_10_2022_clean.R' will draw upon supplied or cited source files to generate and save datasets for analysis in subsequent scripts.

Figure 1:
Figure 1 is a collection of line plots which will help the user visualise how the strength of evolutionary integration between bones within different regions/modules of the avian skeleton change as avian body mass increases. This is achieved by generating many taxonomic subsets of birds with different mean body mass values, and repeating analyses to determin evolutionary integration.

a) Run script 'Residual_Csize_10_08_2022.R' to generate a dataset of avian bone sizes (centroid size computed from landmark constellations), corrected for patterns of allometric scaling (size-dependent scaling; this must be removed because the role of body mass is being investigated and size-dependent scaling could obfuscate such patterns).

b) Run script 'Integration_v_body_mass_11_08_2022.R' to run the analysis.

c) Run script 'Plot_integration_v_body_mass_14_11_2022.R' to plot the analysis output (Figure 1).
