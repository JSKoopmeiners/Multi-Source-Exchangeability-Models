# Multi-Source-Exchangeability-Models

This repository contains the code used to completion the simulation study and data illustration in a manuscript entitled "Dynamic Multi-Resolution Smoothing Using Multi-Source Exchangeability Models" by Kaizer, Koopmeiners and Hobbs. Briefly, multi-source exchangeability models (MEMs) are a flexible Bayesian approach for incorporating supplementary data into the analysis of a primary data source. The following files are included:

MEM_Functions.R: A set of cuntions for implementing MEMs and used for the simulation study and illustration

MEM_Simulations.R: R file to run the simulations discussed in the manuscript for evluating the performance of MEMs

Simulation Results: This is a sub-folder that contains the simulation output from MEM_Simulations.R and used to create the tables and figures in the manuscript.

Figure3_and_SupplementaryFigures_code.R: Code to create Figure 3 in the manuscript and the Figures in the supplementary materials

Figure4_and_MiscellaneousSummaries_code.R: Code to create Figure 4 and miscellaneous summaries reported in the manuscript

Table1_CENICp1Example_code.R: Code for completing real-data example reported in manuscript. Note - the raw data for CENICp1 are not publicaly available but only the sample mean and standard deviation are required to complete the analysis and they are provided in this file

re_model_MetaAnalysis.txt: BRugs code for implementing standard hierarchical model

inits_MetaAnalysis.txt: initial values for standard hierarchical model with three supplementary sources used in simulation and example

inits_MetaAnalysis_2source.txt: initial values for standard hierarchical model with two supplementary sources used in simulation and example
