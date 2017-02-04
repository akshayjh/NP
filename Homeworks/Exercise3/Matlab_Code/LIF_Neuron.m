%% Exercise 3 - Leaky Integrate and Fire Neuron
function [V t] = LIF_Neuron(I, T_max, sinus)


% Neuron parameters:
C_m = 1e-6; % Membrane capacitance
g_leak = 100e-6; % Leakage conductance
V_rest = -60e-3; % Resting potential
V_thr = -20e-3; % Threshold potential
V_spike = 20e-3; % Spike potential

t_0 = 0; % starting time
dt = 25e-6; % step length
t = t_0:dt:T_max;
V = zeros(1,length(t));
V(1) = V_rest; % resting potential as starting point

if sinus == 0
    I_input = zeros(1,length(t));
    I_input(1:(50e-3/dt)) = I;
end
%% rectified sine wave as input I_input(t)
if sinus == 1
A1 = I; % amplitude for rectified sine wave and 10 mikro amps amplitude
f = 50; % frequency for sine wave
I_input = abs(A1 * sin(2*pi*f*t));
I_input(50e-3/dt+1:end) = 0;
%I_input(I_input > 0) = A1;
%I_input(I_input < 0) = 0;
end
%% Calculate LIF model

for i=1:length(t)-1
    if V(i) < V_thr
        V(i+1) = V(i) + dt/C_m * (-g_leak*(V(i) - V_rest) + I_input(i));
    end
    if V(i) >= V_thr
        V(i+1) = V_spike;
    end
    if V(i) == V_spike
        V(i+1) = V_rest;
    end
    
end

% plot model
figure;
plot(t*1000,V*1000, 'r','LineWidth', 1.5, 'Color', 'r'), title('LIF Neuron'), ylabel('V(t) in mV'), xlabel('t in ms');
grid on;

end


