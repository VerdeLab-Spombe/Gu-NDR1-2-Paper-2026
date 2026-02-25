# Gu-NDR1/2-Cdc42-Pard3-Paper-2026

This repository contains code used in the study:
"NDR1/2 kinases regulate cell polarization and cell motility through Cdc42 GTPase and Pard3 signaling in mammalian cells," 
written by Jun Gu, Jelena Marjanovic, Marjana Tomic-Canic, Fangliang Zhang, and Fulvia Verde. Additional information can be found in the Methods section of the paper.

# Contents:

## Angle analysis of migrating cells
The cell X and Y coordinates were obtained using TrackMate in FIJI software. The coordinate information was then analyzed using R scripts provided in 'Angle for migration cells-Fig 2-E-F-I-J'.

## Angle analysis of Golgi orientation between nucleus and Golgi
The X and Y coordinates of the centroids of the nucleus and Golgi marker GM130 were determined using FIJI software. The R scripts in the file 'Golgi orientation-Fig-3-Fig S4' were then used to analyze the angle between the nucleus and Golgi.

## Mean Square Displacement (MSD) analysis
Mean square displacement (MSD) was analyzed using the Sojourner[1] R package according to the manual by uploading ImageJ Particle Tracker files. The related file is 'MSD-with testing samples-Fig 2-K-L-Fig S3 C D'.

[1] Liu S, Yoo S, Tang X, Sung Y, Wu C (2020). Sojourner: statistical analysis of single molecule trajectories. R package version 1.3.0. https://github.com/sheng-liu/sojourner. DOI: 10.18129/B9.bioc.sojourner.

## Cell migrating trajectories
MATLAB scripts in the folder 'Migration-random/wound cell trace' were used for visualization of cell migrating trajectories. The related files are 'Migration-random-Cell trace-Fig S3-A' and 'Migration-wound-Cell trace-Fig2-B'.

## Analysis of the Cdc42 biosensor
MATLAB and R scripts were used for Figures 1 and S2. R and Python scripts were used for data wrangling, while MATLAB scripts were used for analysis of the cross-correlation coefficient between the biosensor and edge velocity, as well as autocorrelation analysis of the periodicity of biosensor activity and edge movements.

For Figures 1 and S2 (Cdc42 biosensor analysis), we first used the u-register software from the Danuser Lab (https://github.com/DanuserLab/u-register). All MATLAB code was run using MATLAB version R2024a.

The related files are:
- 'MATLAB codes-Fig 1 D-J-Fig S2'
- 'Python code_txt_wrangling-Fig 1 D-J-Fig S2'
- 'R code_needed for data analysis-Fig 1 D-J-Fig S2'

Example .csv files are included in each folder.



