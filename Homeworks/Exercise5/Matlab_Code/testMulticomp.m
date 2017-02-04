%% Multicompartement Tests
%% create Stimulus i_stimulus:
clear, clc

t_0 = 0; % starting time in ms
T_max = 100; % max time in ms
dt = 25e-3; % step length in ms
t = t_0:dt:T_max;
%% case 1: create i_stimulus of 5ms and 10 microAmpere/cm2 amplitude
V = [];
i_stimulus = zeros(100, length(t));
i_stimulus(1, 1:length(0:dt:5)) = 10;

V = MultiCompartement(i_stimulus);
%% plot 1
h= pcolor(V);
colormap gray
set(h, 'edgecolor', 'non');
% change axis labels
title('I_{stimulation}: 5ms long pulse, A = 10 \muA/cm^2')
xlabel('t in ms')
ylabel('Compartment #')

%% case 2: 
V = [];
i_stimulus = zeros(100, length(t));
i_stimulus(20, 1:length(0:dt:5)) = 10;
i_stimulus(80, 1:length(0:dt:5)) = 10;

% case 1: test Multicompartement model
V = MultiCompartement(i_stimulus);
%% plot 2
h= pcolor(V);
colormap gray
set(h, 'edgecolor', 'non');
% change axis labels
title('I_{stimulation}: 5ms long pulse, A = 10 \muA/cm^2')
xlabel('t in ms')
ylabel('Compartment #')

%% 
