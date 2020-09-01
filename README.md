# Supplemental Materials to "Algorithmic Jingle Jungle: A Comparison of Implementations of Principal Axis Factoring and Promax Rotation in R and SPSS"

This repository contains the data (i.e., the files we are allowed to share publicly), analysis scripts and figures for the three studies reported in "Algorithmic Jingle Jungle: A Comparison of Implementations of Principal Axis Factoring and Promax Rotation in R and SPSS."


*best_implementations.xlsx* contains a table describing the implementations that performed best in the simulation analyses, as well as indicators of their performances regarding goodness of fit, proportion of occurrence of Heywood cases, and indicator-to-factor correspondences.

*data* contains two subfolders, *HCA*, and *other_datasets*, which contain data files from the Human Cognitive Abilities book by Carroll, and from various different sources (on personality, risk taking etc.), respectively. 

*output* contains four subfolders (*real_data*, *rep_orig*, *simulation_psych_spss*, and *simulation_settings*) with objects created in the analyses. One such object is the *model_recovery_results.RDS* file, which is stored on the top level of the osf.io/a836q page, as it is too large to share via git. To run the R psych vs. SPSS simulation analyses, make sure to download this file separately and place it in *output/simulation_psych_spss*.

*plots* contains pdfs of the different visualizations of the analyses.

*r* contains the R scripts used to run the simulations, analyses, and for generating the plots. Where necessary, the script names are prefixed with a number to indicate the order in which they have to be executed, if the whole analysis sequence should be run. However, as the different output files are available in *output*, the scripts can also be executed without rerunning all previous simulation scripts.

*Real_data_description.xlsx* contains a table with detailed information regarding specifications and sources of the real data sets used in the analyses.