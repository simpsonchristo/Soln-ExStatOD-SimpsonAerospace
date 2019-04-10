clear;
clc;
%Measurement Modeling Example
%Consider a satellite in an equatorial posigrade circular orbit with an
%altitude of 600 km above a spherical Earth. Assume the satellite is 20
%deg. past the zenith direction of a two-way ranging station, which places
%the satellite at 4.3 deg elevation w/r to the station. Assume a signal is
%transmitted from the station at t=0 sec.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simpson Aerospace (c) 2019
%Range rate - Transmitter Beacon

mu = 3.9860044e+14; %m^3/s^2, Earth gravitational parameter
we = (2*pi()/86164);%rad/sec, Earth avg rotational rate
re = 6378137;%meters, spherical Earth radius
alt= 600000;%meters, altitude of satellite
tT1= 0;%sec, time of transmit
dt = 1e-3;%sec, transmit time interval
fT = 24.25e9;%Hz, transmit frequency
c  = 299792458;%m/s, speed of light

%gs vector
phi = acos(cos(deg2rad(20))/cos(deg2rad(4.3)));%rad, latitude of gs
lat = phi;
lon = 0;%rad, longitude of gs

rgsvec  = [re*cos(lat)*1;...
           re*cos(lat)*0;...
           re*sin(lat)];

%initial orbital elements
a    = re+alt;%meters, semimajor axis
ecc  = 0;
inc  = 0;%deg, inclination
raan = 0;%deg, raan, technically NaN
w    = 0;%deg, argument of perigee
nu   = 20;%deg, true anomaly

p = a*(1-ecc^2);%meters, semilatus rectum
h = sqrt(mu*p);%m^2/s, specific angular momentum
E = acos((a/a)*cos(deg2rad(nu)) + ecc);%rad, eccentric anomaly
M = E - ecc*sin(E);%rad, mean anomaly
T = 2*pi()*sqrt(a^3/mu);%sec, period of orbit
n = sqrt(mu/a^3);%rad/sec, mean motion

%pack orbital elements
oe(1) = a;
oe(2) = ecc;
oe(3) = inc;
oe(4) = raan;
oe(5) = w;
oe(6) = rad2deg(M);

%CALCULATE RHO1 AND RHO2
%rho1 = sat at t=tT1 and gs at t=tR1
%rho2 = sat at t=tT2=tT1+dt and gs at t=tR2
[rho1, rho2] = rangefordoppler(oe, nu, tT1, dt, rgsvec);

%fR
N12_dt = (fT/c)*((rho2-rho1)/dt);
fR = fT - N12_dt;

fprintf('Doppler count: %f \n', N12_dt*dt);
fprintf('Received frequency: %f GHz\n', fR/1e9);
if((fT-fR)<fT)
    fprintf('Apparent frequency is lower than actual frequency.\n');
    fprintf('The satellite is moving away from the ground station.\n');
    fprintf('(f_T - f_R)<f_T, (f_T-f_R) = %f GHz\n', (fT-fR)/1e9);
elseif((fT-fR)>fT)
    fprintf('Apparent frequency is greater than actual frequency.\n');
    fprintf('The satellite is moving towards the ground station.\n');
    fprintf('(f_T - f_R)>f_T, (f_T-f_R) = %f GHz\n', (fT-fR)/1e9);
end
rho2 = rho1 + ((c/fT)*N12_dt*dt);
fprintf('Range estimate at t_{R2}: %f', rho2);