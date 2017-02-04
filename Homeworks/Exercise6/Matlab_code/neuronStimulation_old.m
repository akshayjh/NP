%% basic model parameters :
clear;
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
Temp = 6.3; % for debugging purposes
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

%% Vector allocations
T_max = 30;
T_min = 0;
dt = 25e-3;
t_steps = T_min:dt:T_max;
c = 1; % in �F/cm^2

%% Initialisieren von V_m-Vektor
V = zeros(100, length(t_steps));
V(:,1) = V_rest;

V_e = zeros(100,length(t_steps));

m = x_inf(alpha_m(V(:,1)), beta_m(V(:,1)));
n= x_inf(alpha_n(V(:,1)), beta_n(V(:,1)));
h = x_inf(alpha_h(V(:,1)), beta_h(V(:,1)));


%% modified LIF-Neuron:
% define function handles
A_m = @(V,t) -1./tau(alpha_m(V), beta_m(V), k_comp);
B_m = @(V,t) x_inf(alpha_m(V), beta_m(V))./tau(alpha_m(V), beta_m(V), k_comp);
A_n = @(V,t) -1./tau(alpha_n(V), beta_n(V), k_comp);
B_n = @(V,t) x_inf(alpha_n(V), beta_n(V))./tau(alpha_n(V), beta_n(V), k_comp);
A_h = @(V,t) -1./tau(alpha_h(V), beta_h(V), k_comp);
B_h = @(V,t) x_inf(alpha_h(V), beta_h(V))./tau(alpha_h(V), beta_h(V), k_comp);

%% Parameter
rho_axon = 0.009;
r_axon = 1.5e-4;
l_comp = 0.5e-4;

%% Berechnung von R_a
A_axon = pi*r_axon^2;
R_a = rho_axon * (l_comp / A_axon);


%% Initialisieren von C mit Randbedingung Ende abgeschnitten
C = diag(-2*ones(100,1), 0) + diag(ones(99,1), 1) + diag(ones(99,1), -1);
C(1) =-1;
C(end) = -1;
I = eye(100);

%% Felder berechnen
phi = @(rho, x, d, I) rho*I/(4*pi*sqrt(x^2+d^2));
%AF = @(rho, x, d, I) = I*rho*(2*x^2-d^2)/(4*pi*(x^2+d^2)^2.5);

rho = 300;
t_phase = 1;
d = 10e-4;
I_elektrode = zeros(1,length(t_steps));
t_stim = 0:dt:t_phase;
I_elektrode((5/dt):(length(t_stim)+ 5/dt))=-2; % in mA
% plot(I_elektrode) % for test

% biphasic:
I_elektrode((length(t_stim)+5/dt+1):(length(t_stim)+5/dt)+1/dt)=2; % in mA
plot(I_elektrode, 'LineWidth', 2)



%% Berechnung von A_axon = const.
A_axon = I-dt/(c*R_a)*C;




for i=1:length(t_steps)-1
    %% Stimulations-Potential berechnen
    k=1;
    for x=-25e-4:0.5e-4:25e-4
        V_e (k,i) = phi(rho, x, d, I_elektrode(i));
        k=k+1;
    end
    

    %% solve ODEs for m with exponential Euler
    
    m = m .* exp(A_m(V(:,i)) .* dt) + (B_m(V(:,i))/A_m(V(:,i))) * (exp(A_m(V(:,i)) .* dt) - 1);
    n = n .* exp(A_n(V(:,i)) .* dt) + (B_n(V(:,i))/A_n(V(:,i))) * (exp(A_n(V(:,i)) .* dt) - 1);
    h = h .* exp(A_h(V(:,i)) .* dt) + (B_h(V(:,i))/A_h(V(:,i))) * (exp(A_h(V(:,i)) .* dt) - 1);
    
    %% compute ionic current
    i_L = g_L .* (V(:,i) - V_L); % leak current of the cell
    i_K = g_K .* n.^4 .* (V(:,i) - V_K); % potassium ion current
    i_Na = g_Na .* m.^3 .* h .* (V(:,i) - V_Na); % sodium ion current
    i_ion = i_Na + i_K  + i_L; % total ion current
    
    
    %% compute new voltage value with forward Euler
    V(:,i+1) = A_axon \ (V(:,i) + dt/c * (-i_ion + (1/R_a).*C*V_e(1:end-1,i)));
    %V(i+1) = V(i) - 1/C .* (i_ion(i+1)).* dt;
    
    
    
    
end
figure(1)
h=pcolor(V);
% colormap(gray())
c = colorbar
title(['Axon stimulation, mono-phasic, I=5mA'])
set(h,'edgecolor','none');
xlabel('time in \itms')
ylabel('Compartment #')
ylabel(c, 'Action potential in mV')
set(gca,'FontSize',20);

grid on 
% figure(2)
% for i = 1:100
%    plot(V(i,:) + 10*i, 'k')
%    hold on;
% end







