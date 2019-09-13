%% housekeeping

clear
clc
close all

%% info:



% this code is running analysis on motion analysis obtained from SolidWorks
% software for Stirling engine.


%% read data file


% from SolidWorks
FoamDisp = xlsread('Data/Big Lin Disp at CG.xlsx'); % displacement of foam inside controlled volume
HoleAngularDisp = xlsread('Data/Hole angular Displacement.xlsx'); % angular displacement of the hole on the fly wheel 
PistonDisp = xlsread('Data/Small Bottom Face Disp.xlsx'); %displacement of Bottom face of the piston.


% form Vi: engine was tested with 3 different temprature differences 8 10
% 12

EightDeg = load('8degrees_engine3');
TenDeg = load('10degrees_engine3');
TweleveDeg = load('12degrees_engine3');



%% RPM : from angular displacement.

% Get RPM from hole angular displacement.

Angle = HoleAngularDisp(:,3);
Angle = Angle*(pi/180); % convert from deg to radian

%RPM = ( Angle ./ (HoleAngularDisp(:,2)/60)) ./ (2*pi) ;


%% RPM : From linear displacement

