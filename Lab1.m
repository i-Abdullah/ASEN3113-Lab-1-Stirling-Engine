%% housekeeping

clear
clc
close all

%% define constants:

RCylinder = 72 ; % in mm
RFoam = 70 ; % in mm
RPiston = 7.5 ; % in mm
HCylinder = 21; %mm % height of cylinder
HFoam = 11; % mm , height of foam.


%% info:



% this code is running analysis on motion analysis obtained from SolidWorks
% software for Stirling engine.


% add data file to path
addpath('./Data');

%% read data file


% from SolidWorks
FoamDisp = xlsread('Data/Big Lin Disp at CG.xlsx'); % displacement of foam inside controlled volume
HoleAngularDisp = xlsread('Data/Hole angular Displacement.xlsx'); % angular displacement of the hole on the fly wheel 
PistonDisp = xlsread('Data/Small Bottom Face Disp.xlsx'); %displacement of Bottom face of the piston.


% form Vi: engine was tested with 3 different temprature differences 8 10
% 12

T8 = load('8degrees_engine3');
T10 = load('10degrees_engine3');
T12 = load('12degrees_engine3');


%% RPM : From linear displacement

% from linear displacement, the time taken between each two peaks is the
% time needed.

plot(FoamDisp(:,2),FoamDisp(:,3));
[x y] = ginput(2);
period = x(2) - x(1);
frequency = 1/period;
RPM_CAD = frequency * 60 ;

% this RPM matches what we have inputted for solidworks, which is 100;

%% RPM from pressure change
% this's for arctual text.


% RPM for when temp difference is 8
plot(T8(:,1),T8(:,2));
[x y] = ginput(2);
period = x(2) - x(1);
frequency = 1/period;
RPM_T8_Pressure = frequency * 60 ;

% RPM for when temp difference is 10
plot(T10(:,1),T10(:,2));
[x y] = ginput(2);
period = x(2) - x(1);
frequency = 1/period;
RPM_T10_Pressure = frequency * 60 ;

% RPM for when temp difference is 12
plot(T12(:,1),T12(:,2));
[x y] = ginput(2);
period = x(2) - x(1);
frequency = 1/period;
RPM_T12_Pressure = frequency * 60 ;


%% RPM from Sensors

T8_Pass = find(T8(:,8)==1); % see when the wheel passes at the optic, it'll be 1
T10_Pass = find(T10(:,8)==1); % see when the wheel passes at the optic, it'll be 1
T12_Pass = find(T12(:,8)==1); % see when the wheel passes at the optic, it'll be 1

find_T8 = find(diff(T8_Pass)>1,2); % see when the next cycle begins
find_T10 = find(diff(T10_Pass)>1,2); % see when the next cycle begins
find_T12 = find(diff(T12_Pass)>1,2); % see when the next cycle begins

find_T8 = find_T8 + 1; % the diff function reduces index by 1, add that back
find_T10 = find_T10 + 1;
find_T12 = find_T12 + 1;


period = T8(T8_Pass(find_T8(2),1)) - T8(T8_Pass(find_T8(1),1));
RPM_T8_Sensor = (1/period) * 60;

period = T10(T10_Pass(find_T10(2),1)) - T10(T10_Pass(find_T10(1),1));
RPM_T10_Sensor = (1/period) * 60;


period = T12(T12_Pass(find_T12(2),1)) - T12(T12_Pass(find_T12(1),1));
RPM_T12_Sensor = (1/period) * 60;

%% Work: from energy conservation

%% change in volume

Cylinder_Volume = pi*(RCylinder*10^-3)^2 * HCylinder;
Foam_Volume = pi*(RFoam*10^-3)^2 * HFoam;

