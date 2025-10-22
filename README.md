# Hy_CP_Sorting
Electrophysiological analysis and contraction pulse (CP) sorting from MEA recordings in Hydra

## CP Sorting and Feature extraction 

Contraction Pulses (CP) detection and feature extraction were performed from MEA recordings. First and second derivatives were calculated using the centered finite difference method for equally spaced data.
<!--![alt tag](images/HY_Rec.png)-->

## CP Sorting 

Temporal distribution (Raster Plot) of detected peaks in the original full-length traces 
![alt tag](images/RasterPlot_1.png)

## Feature extraction 

Hy CP Sorting computes:
• the single CP duration (CPI), their sum for burst (BurstTime) and their sum
for recording file (CTime);
• the intercontraction burst duration (IcBI) and their sum (ETime);
• the number of bursts (nBurst) and the number of IcBI (nIcBI);
• the number of CP events for each burst (nCPBurst) and recording file (nCP).
![alt tag](images/CPI.png)
![alt tag](images/IcBI.png)
![alt tag](images/Sum_CPburst.png)
![alt tag](images/N_CPburst.png)

## Comparative analysis of data

Traces alignment to the duration of the shortest one to ensure consistent comparisons between recordings.

![alt tag](images/Overlap.png)

## Raster plot

Temporal distribution (Raster Plot) of detected peaks in the original full-length traces (black traces) respect to their analyzed segments (Dotted lines represent the time of the shorter traces)
![alt tag](images/RasterPlot_2.png)


