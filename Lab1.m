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
Cylinder_Volume = pi*(RCylinder*10^-3)^2 * HCylinder*10^-3;
Foam_Volume = pi*(RFoam*10^-3)^2 * HFoam*10^-3;
R_air = 0.287; % (KJ / Kg-K)
Cp_air = 1.005 ; % KJ/(Kg*K)
Cv_air = 0.718 ; % KJ/(Kg*K)
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


V1 = Cylinder_Volume - Foam_Volume;

% we will use piston displacement to compute the change in volume
% we will thus zero the initial displacement,

% the zero point is when the piston is at the bottome, so we will make the
% max position equatl to 0.

PistonDisp_callibrated = PistonDisp(:,3) - min(PistonDisp(:,3)) ; 
% min here will give the maximum in negartive.

% Change in volume
DV = (PistonDisp_callibrated)*10^-3 * (pi*(RPiston*10^-3)^2);

% now to get the volume at any given point by adding V1 to DV.

% one cycle is sufficient,

V2 = max(DV) + V1 ;

%% work: info


% 1 -> 2 : Isobaric Heating
% 2 -> 3 : Isothermal expansion
% 3 -> 4 : Isobaric cooling
% 4 -> 1 : Isothermal compression


%% Work: idealized case

% there are 4 processes


% 4th column is bottom of top face, and 5th is top of bottom, we're taking
% the average temp.

% note: *6.89476 is to convert from psi to kpa.

% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

m_air = ( max((T8(:,2)*6.89476)) .* min((V1+DV)) ) ./ ( R_air .* max(( T8(:,4) + T8(:,5))./2)+273.15); % PV / RT


% 1 -> 2 :

W12_8 = 0;
Qnet_12_8 = m_air*Cv_air*((8+273)) ;

% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

% 2 -> 3 :

W23_8 = m_air*R_air*(max(( T8(:,4) + T8(:,5))./2)+273.15)*log(V2/V1) ;
Qnet_23_8 = W23_8 ;

% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

% 3 -> 4 :

W34_8 = 0 ;
Qnet_34_8 = -m_air*Cv_air*((8+273)) ;

% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

% 4 -> 1 :

W41_8 =  m_air*R_air*(min(( T8(:,4) + T8(:,5))./2)+273.15)*log(V1/V2);
Qnet_41_8 = -1*m_air*Cv_air*((8+273)) ;

