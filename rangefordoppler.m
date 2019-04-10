function [rho1, rho2] = rangefordoppler(oe,nu,tT1,dt,rgs)
%RANGEFORDOPPLER Determine range between sat and gs
%Range is determined between the satellite position at tT# and the ground
%station position at tR#. Position is determined using a spherical Earth
%and two-body approximation. 
%INPUTS:
%       oe, orbital elements [a; ecc; inc; raan; aop; M0]
%       nu, true anomaly
%       tT1, first transmit time
%       dt, interval between transmissions from beacon
%       rgs, ecf postion of ground station
%OUTPUTS:
%       rho1, range between sat position at tT1 and gs position at tR1
%       rho2, range between sat position at tT2 and gs position at tR2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simpson Aerospace (c) 2019
%Christopher R. Simpson

mu = 3.9860044e+14; %m^3/s^2, Earth gravitational parameter
we = (2*pi()/86164);%rad/sec, Earth avg rotational rate
c  = 299792458;%m/s, speed of light

%RHO1
%xt = ecf position of sat at tT1
XtT1 = oe2rv(oe,nu);
agT1 = we*tT1;%rad, Greenwich angle at tT1
rotz = [cos(agT1) -sin(agT1) 0; sin(agT1) cos(agT1) 0; 0 0 1];
xtT1 = XtT1(1,1:3)*rotz;
instantrho = norm(xtT1 - rgs);

tR1 = instantrho/c;
agR1 = we*tR1 + agT1;%rad, Greenwich angle at tR1
rotz = [cos(agR1) -sin(agR1) 0; sin(agR1) cos(agR1) 0; 0 0 1];
rho1 = norm(XtT1(1,1:3) - (transpose(rotz)*rgs));

%RHO2
tT2  = dt + tT1;%sec, second transmit time
a    = oe(1);%meters, semimajor axis
ecc  = oe(2);%   , eccentricity

EtT1 = acos((a/a)*cos(deg2rad(nu)) + ecc);%rad, eccentric anomaly
MtT1 = EtT1 - ecc*sin(EtT1);%rad, mean anomaly
n = sqrt(mu/a^3);%rad/sec, mean motion

EtT2 = keplerseqn(ecc,n,dt,MtT1);%rad, eccentric anomaly at ta
nutT2 = acos((cos(EtT2) - ecc)/(1 - ecc*cos(EtT2)));%rad, true anomaly at ta
oe(6) = n*tT2 + MtT1;

XtT2 = oe2rv(oe,nutT2);
agT2 = we*tT2 + agT1;%rad, Greenwich angle at tT1
rotz = [cos(agT2) -sin(agT2) 0; sin(agT2) cos(agT2) 0; 0 0 1];
xtT2 = XtT2(1,1:3)*rotz;
instantrhoT2 = norm(xtT2 - rgs);

tR2 = instantrhoT2/c;
agR2 = we*tR2 + agT2;%rad, Greenwich angle at tR1
rotz = [cos(agR2) -sin(agR2) 0; sin(agR2) cos(agR2) 0; 0 0 1];
rho2 = norm(XtT2(1,1:3) - (transpose(rotz)*rgs));

end
