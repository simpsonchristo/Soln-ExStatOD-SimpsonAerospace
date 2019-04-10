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
%Two-way range estimate

mu = 3.9860044e+14; %m^3/s^2, Earth gravitational parameter
we = (2*pi()/86164);%rad/sec, Earth avg rotational rate
re = 6378137;%meters, spherical Earth radius
alt= 600000;%meters, altitude of satellite
t_t= 0;%sec, time of transmit
c  = 299792458;%m/s, speed of light

%INSTANTANEOUS RANGE
%1 - dot(zthat, rvec) = cos(deg2rad(4.3));
%2 - cos(phi)*(dot(zthat,rvec)) = cos(deg2rad(20));
phi = acos(cos(deg2rad(20))/cos(deg2rad(4.3)));%rad, latitude of gs
lat = phi;
lon = 0;%rad, longitude of gs

rsatvec = [(re+alt);...
           0;...
           0];
rgsvec  = [re*cos(lat)*1;...
           re*cos(lat)*0;...
           re*sin(lat)];
rho = rsatvec - rgsvec;%meters, range differenced in ecf

%SIGNAL ARRIVAL TIME
ta = t_t + norm(rho)/c;%sec, estimated time of arrival

%NEW RANGE
%determine orbital elements
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

Ma = n*(ta);%rad, mean anomaly at ta
Ea = keplerseqn(ecc,n,ta,M);%rad, eccentric anomaly at ta
nu_a = acos((cos(Ea) - ecc)/(1 - ecc*cos(Ea)));%rad, true anomaly at ta
%pack orbital elements
oe_a(1) = a;
oe_a(2) = ecc;
oe_a(3) = inc;
oe_a(4) = raan;
oe_a(5) = w;
oe_a(6) = rad2deg(Ma);

Xt = oe2rv(oe_a,nu_a);
ag = we*ta;
rotz = [cos(ag) -sin(ag) 0; sin(ag) cos(ag) 0; 0 0 1];
xt = Xt(1,1:3)*rotz;
xt = transpose(xt);

rho_a = xt - rgsvec;

fprintf('Range at t=0: %f meters\n', norm(rho));
fprintf('Range at time of arrival, t = %f: %f meters\n',ta, norm(rho_a));




