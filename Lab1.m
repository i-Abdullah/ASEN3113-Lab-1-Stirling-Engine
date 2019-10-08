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
gamma_air = 1.4 ; % cv/cp
Cp_air = 1.005 ; % KJ/(Kg*K)
Cv_air = 0.718 ; % KJ/(Kg*K)
Cv_Air = R_air / (gamma_air - 1) ;


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


%FoamDisp_ = xlsread('Data/RPM_Data/Big Lin Disp at CG.xlsx'); % displacement of foam inside controlled volume
%HoleAngularDisp = xlsread('Data/RPM_Data/Hole angular Displacement.xlsx'); % angular displacement of the hole on the fly wheel 
PistonDisp_8 = xlsread('Data/RPM_Data/SmallBottom96.xlsx'); %displacement of Bottom face of the piston.
PistonDisp_10 = xlsread('Data/RPM_Data/SmallBottom112.xlsx'); %displacement of Bottom face of the piston.
PistonDisp_12 = xlsread('Data/RPM_Data/SmallBottom129.xlsx'); %displacement of Bottom face of the piston.

% zero data:

PistonDisp_8(:,3) = PistonDisp_8(:,3) - min(PistonDisp_8(:,3));
PistonDisp_10(:,3) = PistonDisp_10(:,3) - min(PistonDisp_10(:,3));
PistonDisp_12(:,3) = PistonDisp_12(:,3) - min(PistonDisp_12(:,3));

% make sure that we start from displacement = 0:

j = find(PistonDisp_8(:,3)==min(PistonDisp_8(:,3))); % 
PistonDisp_8(1:j-1,:) = [];

j = find(PistonDisp_10(:,3)==min(PistonDisp_10(:,3))); % 
PistonDisp_10(1:j-1,:) = [];

j = find(PistonDisp_12(:,3)==min(PistonDisp_12(:,3))); % 
PistonDisp_12(1:j-1,:) = [];


% zero time

PistonDisp_8(:,2) = PistonDisp_8(:,2) - (PistonDisp_8(1,2));

PistonDisp_10(:,2) = PistonDisp_10(:,2) - (PistonDisp_10(1,2));

PistonDisp_12(:,2) = PistonDisp_12(:,2) - (PistonDisp_12(1,2));


% form Vi: engine was tested with 3 different temprature differences
% 8 10 12

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
% this's for actual experiment.


%{


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


%}

%% RPM from Sensors

% there's optical sensor that checks when the 
T8_Pass = find(T8(:,8)==1); % see when the wheel passes at the optic, it'll be 1
T10_Pass = find(T10(:,8)==1); % see when the wheel passes at the optic, it'll be 1
T12_Pass = find(T12(:,8)==1); % see when the wheel passes at the optic, it'll be 1

find_T8 = find(diff(T8_Pass)>1,2); % see when the next cycle begins
find_T10 = find(diff(T10_Pass)>1,2); % see when the next cycle begins
find_T12 = find(diff(T12_Pass)>1,2); % see when the next cycle begins

find_T8 = find_T8 + 1; % the diff function reduces index by 1, add that back
find_T10 = find_T10 + 1;
find_T12 = find_T12 + 1;



% those fin_T8, T10, T10 are the two indices, which are the indices that
% represent the beggining and end of the cycle. from them we can get
% period, and information about all the cyclces.

period = T8(T8_Pass(find_T8(2)),1) - T8(T8_Pass(find_T8(1)),1);
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

%% calibrating two data's

% the issue is that experimntal data isn't the same frequency as the actual
% data.

% we can figure out how much time it takes for one cycle from the rpm data,
% this should be the same for both CAD model and Experiment.

% get time it took for one cycle:


% col 1 is time:
Start_t_8 = 0;
End_t_8 = T8(T8_Pass(find_T8(2),1)) - T8(T8_Pass(find_T8(1),1)); % Time of each cycle is constant so we can assume we start from 0 and go to 0 + dt.

Start_t_10 = 0;
End_t_10 = T10(T10_Pass(find_T10(2),1)) - T10(T10_Pass(find_T10(1),1));

Start_t_12 = 0;
End_t_12 = T12(T12_Pass(find_T12(2),1)) - T12(T12_Pass(find_T12(1),1));


% pull out all data realted to this time period between the start and the
% end, it'll not be the same length, so we will need to pull that to adjust
% it

i = knnsearch(PistonDisp_8(:,2),End_t_8); % find where times are the same.

Volume_8 = V1 + PistonDisp_8(1:i,3)*10^-3 * (pi*(RPiston*10^-3)^2) ;
t_Volume_8 = PistonDisp_8(1:i,2);


i = knnsearch(PistonDisp_10(:,2),End_t_10); % find where times are the same.
i = i + 200 ;
Volume_10 = V1 + PistonDisp_10(1:i,3)*10^-3 * (pi*(RPiston*10^-3)^2) ;
t_Volume_10 = PistonDisp_10(1:i,2);


i = knnsearch(PistonDisp_12(:,2),End_t_12); % find where times are the same.
i = i + 200 ;
Volume_12 = V1 + PistonDisp_12(1:i,3)*10^-3 * (pi*(RPiston*10^-3)^2) ;
t_Volume_12 = PistonDisp_12(1:i,2);



Pressure_8 = T8(T8_Pass(find_T8(1)):T8_Pass(find_T8(2)),2)  ; 
Pressure_8 = Pressure_8 + 12 ; % add atmospheric pressur in psi
Pressure_8 = Pressure_8* 6894.76 ; % convert to Pa

Pressure_10 = T8(T8_Pass(find_T8(1)):T8_Pass(find_T8(2)),2)  ; 
Pressure_10 = Pressure_10 + 12 ; % add atmospheric pressur in psi
Pressure_10 = Pressure_10* 6894.76 ; % convert to Pa


Pressure_12 = T8(T8_Pass(find_T8(1)):T8_Pass(find_T8(2)),2)  ; 
Pressure_12 = Pressure_12 + 12 ; % add atmospheric pressur in psi
Pressure_12 = Pressure_12* 6894.76 ; % convert to Pa


t_Pressure_8 = T8(T8_Pass(find_T8(1)):T8_Pass(find_T8(2)),1) ;
t_Pressure_8 = t_Pressure_8 - t_Pressure_8(1) ; 
Volume_8_interp = interp1(t_Volume_8,Volume_8,t_Pressure_8);

t_Pressure_10 = T8(T8_Pass(find_T8(1)):T8_Pass(find_T8(2)),1) ;
t_Pressure_10 = t_Pressure_10 - t_Pressure_10(1) ; 
Volume_10_interp = interp1(t_Volume_10,Volume_10,t_Pressure_10);


t_Pressure_12 = T8(T8_Pass(find_T8(1)):T8_Pass(find_T8(2)),1) ;
t_Pressure_12 = t_Pressure_12 - t_Pressure_12(1) ; 
Volume_12_interp = interp1(t_Volume_12,Volume_12,t_Pressure_12);


figure(1)
scatter(Volume_8_interp,Pressure_8)
grid minor
xlabel('Volume [m^3]');
ylabel('Pressure [Pa]');
title(['PV diagram for T = 8' char(176) ' C stirling engine'])

 
figure(2)
scatter(Volume_10_interp,Pressure_10)
grid minor
xlabel('Volume [m^3]');
ylabel('Pressure [Pa]');
title(['PV diagram for T = 10' char(176) ' C stirling engine'])


figure(3)
scatter(Volume_12_interp,Pressure_12)
grid minor
xlabel('Volume [m^3]');
ylabel('Pressure [Pa]');
title(['PV diagram for T = 12' char(176) ' C stirling engine'])

hi = 0;




% now make them the same size:






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
m_air = 9;

% 1 -> 2 :

W12_8 = 0;
Qnet_12_8 = m_air*Cv_air*((8)) ;

W12_10 = 0;
Qnet_12_10 = m_air*Cv_air*((10)) ;

W12_12 = 0;
Qnet_12_12 = m_air*Cv_air*((12)) ;


% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

% 2 -> 3 :

W23_8 = m_air*R_air*(max(( T8(:,4) + T8(:,5))./2)+273.15)*log(V2/V1) ;
Qnet_23_8 = W23_8 ;

W23_10 = m_air*R_air*(max(( T10(:,4) + T10(:,5))./2)+273.15)*log(V2/V1) ;
Qnet_23_10 = W23_10 ;

W23_12 = m_air*R_air*(max(( T12(:,4) + T12(:,5))./2)+273.15)*log(V2/V1) ;
Qnet_23_12 = W23_12 ;


% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

% 3 -> 4 :

W34_8 = 0 ;
Qnet_34_8 = -m_air*Cv_air*((8)) ;

W34_10 = 0 ;
Qnet_34_10 = -m_air*Cv_air*((10)) ;

W34_12 = 0 ;
Qnet_34_12 = -m_air*Cv_air*((12)) ;


% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

% 4 -> 1 :

W41_8 =  m_air*R_air*(min(( T8(:,4) + T8(:,5))./2)+273.15)*log(V1/V2);
Qnet_41_8 = W41_8 ;

W41_10 =  m_air*R_air*(min(( T10(:,4) + T10(:,5))./2)+273.15)*log(V1/V2);
Qnet_41_10 = W41_10 ;

W41_12 =  m_air*R_air*(min(( T12(:,4) + T12(:,5))./2)+273.15)*log(V1/V2);
Qnet_41_12 = W41_12 ;


% themral efficiency: 
% ((W23_8 + W41_8)/(Qnet_12_8+Qnet_23_8))*100


% ideal efficiency from stilring

T1_8 = (min(( T8(:,4) + T8(:,5))./2)+273.15);
T2_8 = T1_8;
T3_8 = (max(( T8(:,4) + T8(:,5))./2)+273.15);
T4_8 = T3_8;

T1_10 = (min(( T10(:,4) + T10(:,5))./2)+273.15);
T2_10 = T1_10;
T3_10 = (max(( T10(:,4) + T10(:,5))./2)+273.15);
T4_10 = T3_10;

T1_12 = (min(( T12(:,4) + T12(:,5))./2)+273.15);
T2_12 = T1_12;
T3_12 = (max(( T12(:,4) + T12(:,5))./2)+273.15);
T4_12 = T3_12;



nth_ideal_8 = ( T3_8 - T1_8 ) / ( T3_8 + (( Cv_Air *( T3_8 - T2_8 ) ) / R_air*log(V2/V1) ) ) * 100 ;
nth_actual_8 = ( trapz(Volume_8_interp(1:1020),Pressure_8(1:1020)) ./ (Qnet_12_8+Qnet_23_8) ) * 100  ;

nth_ideal_10 = ( T3_10 - T1_10 ) / ( T3_10 + (( Cv_Air *( T3_10 - T2_10 ) ) / R_air*log(V2/V1) ) ) * 100 ;
nth_actual_10 = ( trapz(Volume_10_interp(1:1020),Pressure_10(1:1020)) ./ (Qnet_12_10+Qnet_23_10) ) * 100  ;


nth_ideal_12 = ( T3_12 - T1_12 ) / ( T3_12 + (( Cv_Air *( T3_12 - T2_12 ) ) / R_air*log(V2/V1) ) ) * 100 ;
nth_actual_12 = ( trapz(Volume_12_interp(1:1020),Pressure_12(1:1020)) ./ (Qnet_12_12+Qnet_23_12) ) * 100  ;

%% output results:

Actual = { nth_actual_8 ; nth_actual_10 ; nth_actual_12 } ;
Ideal = { nth_ideal_8 ; nth_ideal_10 ; nth_ideal_12 } ;
Results = { 'T = 8';'T = 10';'T = 12'};
table(Results,Actual,Ideal)