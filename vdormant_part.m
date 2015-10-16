function dydt=vdormant_part(t,y,info)
% function dydt=vdormant_part(t,y,info)
%
% Simulates the S-V system

% Variables
S = y(1);
V = y(2);

% Dynamics
dSdt = -info.phi*S*V*(1+info.delta);
dVdt = -info.phi*S*V;

% Return
dydt = [dSdt; dVdt];
