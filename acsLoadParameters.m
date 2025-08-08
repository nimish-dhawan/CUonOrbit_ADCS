%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NIMISH DHAWAN
% acsLoadParameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters and Design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =======================================================================
% CONTROL GAINS
option = 1;         % 0 for aggressive; 1 for realistic
 
if option == 0
    K_d = 1e-03;
    K_p = 1e-04;    % 5e-5 is a bit more aggressive

elseif option == 1
    K_d = 7.5e-03;
    K_p = 1.25e-04;     

elseif option == 2
    K_d = 0.008;
    K_p = 0.003;
end

K_w = 50e4;

% =======================================================================
% SIM CONFIG
tspan               = 3600 * 10;
write_output        = false;                        % Generate .a attitude file for STK
outputName          = 'attitude_SSO_attTable.a';  
simFile             = 'ADCS_onOrbit_14_4_3.slx';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Date and Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YYYY = 2025;
DD = 01;
MM = 01;
HH = 16;
MIN = 00;
SEC = 000;
dYear = decyear(YYYY,MM,DD);
JD0 = juliandate(YYYY,MM,DD,HH,MIN,SEC);

rE = 6378.14;                        % Earth radius, km
wE = [0 0 15.04*(pi/180)*(1/3600)];  % Earth rotation speed, rad/s
mu = 398600.4418;                    % Gravitational parameter, km^3/s^2
p = 0.00000451;                      % Solar pressure constant, N/m^2
CD = 2.3;                            % Drag Coefficent
eta = 0.3;                           % Absorption Coefficient
mag_moment = [0.005; 0.005; 0.005];  % Magnetic Dipole

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spacecraft and Wheels Initial Scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =======================================================================
% INITIAL CONDITIONS
w_b_ini = [0; 0.1; 0.1];
q_ini   = [0; 0; 0; 1];     % Initial quaternion
m       = 4;                % Spacecraft mass in kg

% =======================================================================
% WHEELS SETUP

A_w     = eye(3);                   % Each wheel aligned with an axis of the s/c     
w_w_ini = [0; 0; 0];                % Wheels initial spin rates in rad/s for wheel 1, 2 and 3 respectively
r_w     = [0 0 0; 0 0 0; 0 0 0];    % Location of reaction wheels from the COM

J_wheels = [12894 0 0; 0 12774 0; 0 0 18389] * 1e-09;     % Cubespace Bigger Wheels [CW0057]

% =======================================================================
% MAGNETORQUER SETUP

r_m = [0 0 0; 0 0 0; 0 0 0];        % Location of the magnetorquers from COM
A_m = eye(3);           % Each magnetorquer aligned with an axis of the s/c

% =======================================================================
% SPACECRAFT BUS SETUP
xAxis = 10e-2;     % in m
yAxis = 10e-2;     % in m
zAxis = 30e-2;     % in m

I_x = (1/12) * m * (yAxis^2 + zAxis^2);
I_y = (1/12) * m * (xAxis^2 + zAxis^2);
I_z = (1/12) * m * (xAxis^2 + yAxis^2);

J_body  = [I_x 0 0; 0 I_y 0; 0 0 I_z];                              % Spacecraft inertia matrix in kg m2
n_list  = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];            % Surface normal vectors
A_list  = [10*30; 10*30; 10*10; 10*30; 10*30; 10*10] * 1e-04;       % Defining panel areas in m2
COP = [xAxis/2   0          0;         % +X face
        0        yAxis/2    0;         % +Y face
        0        0          zAxis/2;   % +Z face
       -xAxis/2  0          0;         % -X face
        0       -yAxis/2    0;         % -Y face
        0        0         -zAxis/2];  % -Z face                    % Assuming COP is at the center of each panel in m. 

% =======================================================================
% KALMAN FILTER INITIALIZATION
P_0 = 0.0000001 * eye(7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensor Modelling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_S_b = 0.05 * [1; 1; 1];
sigma_B_b = 0.05 * [1; 1; 1];
n_drift_0 = 0.0001 * [1; 1; 1] * 0;
w_bias_0 = 0.0001 * [1; 1; 1] * 0;
sigma_n = (0.05) * [1; 1; 1] * 0;
sigma_drift = 0.05 * [1; 1; 1] * 0;
tau = 300;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbital Elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_apogee        = 520;                % km
semimajor_axis  = rE + h_apogee;      % km
e               = 6.91944e-15;
inc             = 97.4065;            % deg
raan            = 132.112;            % deg
omega           = 0;                  % deg
true_anomaly    = 1.42146e-14;        % deg

period = sqrt((4* pi^2 * semimajor_axis^3)/(mu));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Determination and Guidance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get_input = 0;      % Put 0 for OREKIT output and 1 for STK output
% 
% if get_input == 1
%     % Getting the attitude profile from STK output
%     q_nadir_data = readtable('C:\Users\nimis\Documents\STK 12\COonOrbit Scenarios\CUsat\csv data\CUsat_SSO_Attitude_Quaternions_nadir.csv');
%     q_nadir = [q_nadir_data.q1 q_nadir_data.q2 q_nadir_data.q3 q_nadir_data.q4];
% 
%     % Getting angular body rates from STK output
%     w_nadir_data = readtable('C:\Users\nimis\Documents\STK 12\COonOrbit Scenarios\CUsat\csv data\CUsat_SSO_Body_Angular_Rates_nadir.csv');
%     w_nadir = [w_nadir_data.x_rad_sec_, w_nadir_data.y_rad_sec_, w_nadir_data.z_rad_sec_];
% else
% 
%     % Getting data from OREKIT output file
%     nadir_data = readtable('C:\Users\nimis\nadir_attitude_quaternions_log.csv');
%     sun_data = readtable('C:\Users\nimis\sun_attitude_quaternions_log.csv');
% 
%     % Attitude Profiles
%     q_nadir = [nadir_data.q1 nadir_data.q2 nadir_data.q3 nadir_data.q4];
%     q_sun = [sun_data.q1 sun_data.q2 sun_data.q3 sun_data.q4];
% 
%     % Angular Velocities
%     w_nadir = [nadir_data.wX nadir_data.wY nadir_data.wZ];
%     w_sun = [sun_data.wX sun_data.wY sun_data.wZ];
% 
%     % Position
%     pos = [nadir_data.X_m_, nadir_data.Y_m_, nadir_data.Z_m_];
% 
%     % Velocity
%     vel = [nadir_data.vX_m_s_, nadir_data.vY_m_s_, nadir_data.vZ_m_s_];
% 
% end
% 
% % Create time vector (1-second steps)
% t = (0:height(q_nadir)-1)'; 
% 
% % Desired Attitude Profile Input
% % Nadir Pointing
% q_input_nadir.time = t;
% q_input_nadir.signals.values = q_nadir; 
% q_input_nadir.signals.dimensions = 4;
% 
% % Sun Pointing
% q_input_sun.time = t;
% q_input_sun.signals.values = q_sun;
% q_input_sun.signals.dimenstions = 4;
% 
% % Desired Body Angular Rates 
% % Nadir Pointing
% w_input_nadir.time = t;
% w_input_nadir.signals.values = w_nadir; 
% w_input_nadir.signals.dimensions = 3;
% 
% % Sun Pointing
% w_input_sun.time = t;
% w_input_sun.signals.values = w_sun;
% w_input_sun.signals.dimensions = 3;
% 
% if get_input == 0
%     % Position Input
%     pos_input.time = t;
%     pos_input.signals.values = pos;
%     pos_input.signals.dimensions = 3;
% 
%     % Velocity Input
%     vel_input.time = t;
%     vel_input.signals.values = vel;
%     vel_input.signals.dimensions = 3;
% end

% clc

disp('Parameters set')

%%
run('acsRunSim.m')