function [i_ion, gate] = hh_current(V, dt, last_gate)
% Output arguments:
% i_ion : ion current
% gate : current gate state
% Input arguments:
% V : current Voltage V
% dt : time step dt
% last_gate : last gate state state

%% model parameters:
V_Na = 115; % in mV
V_K = -12;
V_L = 10.6; 
V_rest = -70;

g_Na = 120; % in mS/cm^2
g_K = 36;
g_L = 0.3;

C = 1; % in µF/cm^2
%% 



i_L = g_L * (V - V_L); % leak current of the cell 
i_K = g_K * n.^4 .* (V - V_K); % potassium ion current
i_Na = g_Na * m.^3 .* h .* (V - V_Na); % sodium ion current 
i_ion = i_Na + i_K + i_L; % total ion current


end
