function [V, t] = HodgkinHuxley_2(i_stimulus, Temp)
%% basic model parameters :
% alpha
alpha_m = @(V) (2.5 -0.1.*V)./(exp(2.5 - 0.1.*V) - 1); 
alpha_n = @(V) (0.1 - 0.01.*V)./(exp(1-0.1.*V) - 1);
alpha_h = @(V) 0.07.*exp(-V/20);
% beta
beta_m = @(V) 4*exp(-V/18);
beta_n = @(V) 0.125 * exp(-V/80);
beta_h = @(V) 1./(exp(3-0.1*V) +1);

% temperature correction factor k(T in degree Celsius)
k = @(T) 3^(0.1*(T-6.3));
% Temp = 6.3; % for debugging purposes
k_comp = k(Temp);

% function handles to tau and x_inf
tau = @(alpha, beta, K) 1./(K .* (alpha + beta));
x_inf = @(alpha, beta) alpha./(alpha + beta);

V_Na = 115; % in mV
V_K = -12;
V_L = 10.6; 
V_rest = 0;

g_Na = 120; % in mS/cm^2
g_K = 36;
g_L = 0.3;

C = 1; % in µF/cm^2

V(1) = V_rest;
m(1) = x_inf(alpha_m(V(1)), beta_m(V(1)));
n(1) = x_inf(alpha_n(V(1)), beta_n(V(1)));
h(1) = x_inf(alpha_h(V(1)), beta_h(V(1)));
%% Vector allocations
t_0 = 0; % starting time in ms
T_max = 100; % max time in ms
dt = 0.01; % step length in ms
t = t_0:dt:T_max;
V = zeros(1,length(t));
V(1) = V_rest; % resting potential as starting point
%define ODE for exp. Euler:
DGL_exp = @(V, t) A(t)*V + B(t);

%% modified LIF-Neuron:
% define function handles
A_m = @(V,t) -1./tau(alpha_m(V), beta_m(V), k_comp);
B_m = @(V,t) x_inf(alpha_m(V), beta_m(V))./tau(alpha_m(V), beta_m(V), k_comp); 
A_n = @(V,t) -1./tau(alpha_n(V), beta_n(V), k_comp);
B_n = @(V,t) x_inf(alpha_n(V), beta_n(V))./tau(alpha_n(V), beta_n(V), k_comp); 
A_h = @(V,t) -1./tau(alpha_h(V), beta_h(V), k_comp);
B_h = @(V,t) x_inf(alpha_h(V), beta_h(V))./tau(alpha_h(V), beta_h(V), k_comp); 

for i=1:length(t)-1
    %% solve ODEs for m with exponential Euler
    
    m(i+1) = m(i) * exp(A_m(V(i)) .* dt) + (B_m(V(i))/A_m(V(i))) * (exp(A_m(V(i)) .* dt) - 1);
    n(i+1) = n(i) * exp(A_n(V(i)) .* dt) + (B_n(V(i))/A_n(V(i))) * (exp(A_n(V(i)) .* dt) - 1);
    h(i+1) = h(i) * exp(A_h(V(i)) .* dt) + (B_h(V(i))/A_h(V(i))) * (exp(A_h(V(i)) .* dt) - 1);
    
    %% compute ionic current
    i_L(i+1) = g_L * (V(i) - V_L); % leak current of the cell 
    i_K(i+1) = g_K * n(i+1)^4 * (V(i) - V_K); % potassium ion current
    i_Na(i+1) = g_Na * m(i+1)^3 .* h(i+1) .* (V(i) - V_Na); % sodium ion current 
    i_ion(i+1) = i_Na(i+1) + i_K(i+1) + i_L(i+1); % total ion current
    
    %% compute new voltage value with forward Euler
    V(i+1) = V(i) - 1/C .* (i_ion(i+1) - i_stimulus(i+1)).* dt;
    %V(i+1) = V(i) - 1/C .* (i_ion(i+1)).* dt;

    
 
    
end

end