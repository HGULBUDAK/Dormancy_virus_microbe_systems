function dydt=vdormant_full(t,y,info)
% function dydt=vdormant_full(t,y,info)
%
% Simulates the S-C-D-I-V system

% Variables
S = y(1);
C = y(2);
D = y(3);
I = y(4);
V = y(5);

% Dynamics
dSdt = -info.kplus*S*V+(1-info.p)*info.kminus*C;
dCdt = info.kplus*S*V-info.kminus*C-info.kforw*C;
dDdt = info.kminus*info.p*C;
dIdt = info.kforw*C;
dVdt = -info.kplus*S*V+info.kminus*C;

% Return
dydt = [dSdt; dCdt; dDdt; dIdt; dVdt];
