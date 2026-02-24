# Gu-NDR1/2-Cdc42-Pard3-Paper-2026
This repository contains code used in NDR1/2 kinases regulate cell polarization and cell motility through Cdc42 GTPase and Pard3 signaling in mammalian cells, written by Jun Gu, Jelena Marjanovic, Marjana Tomic-Canic, Fangliang Zhang, and Fulvia Verde. Additional information of this paper can be found in the methods section. 
# Contents:
## Angel analysis of migrating cells.
The cell X an Y corordinate was obtained via TrackMate in FIJI software and the the information of corordinate was applied in R scripts as shown in Angel for migration cells-Fig 2-E-F-I-J.
## Angel analysis of Golgi orientation between  nuclear and Golgi.
The X, Y corordiante of the centroid of nuclear and Golgi marker GM130 was determined by FIJI software. Then the R scripts in the folder Golgio orientation-Fig-3-Fig S4 was used to analyze the angel between nuclear and GOlgi.
## Mean Square Displacement (MSD) analysis.
The mean square displacement (MSD) was analyzed using Sojourner R package according the manual by upload ImageJ Particle Tracker file.
## Cell migrating trajectories 
The Matlab scripts in the folder Migration-random/wound cell trace were used for the data visulization of cell migrating trajectories.
## The analysis of Cdc42 biosensor
The Matlab scripts and R scripts in Fig 1 and Fig S2 were used. 1. The R and python scripts were used for the data wrangling and Matlab scripts used for analysis of crosscoefficiency between Biosensor and edge velocity, and autocorrelation analysis the periodity of biosensor and edge movements. The Fig 1 and Fig S2 Cdc42 biosensor analysis, we first used the u-registration software from Danuser Lab software page from github. All Matlab code was run using Matlab Version R2024a.



