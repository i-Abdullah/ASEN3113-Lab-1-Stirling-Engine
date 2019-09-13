%% housekeeping

clear
clc
close all

%% info:



% this code is running analysis on motion analysis obtained from SolidWorks
% software for Stirling engine.


%% read data file



FoamDisp = xlsread('/Data/Big Lin Disp at CG.xlsx'); % displacement of foam inside controlled volume
HoleAngularDisp = xlsread('/Data/Hole angular Displacement.xlsx'); % angular displacement of the hole on the fly wheel 
PistonDisp = xlsread('/Data/Small Bottom Face Disp.xlsx'); %displacement of Bottom face of the piston.


