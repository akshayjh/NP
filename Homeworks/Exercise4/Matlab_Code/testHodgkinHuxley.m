%% test file for Hodgkin & Huxley model neuron
t_0 = 0; % starting time in ms
T_max = 100; % max time in ms
dt = 0.01; % step length in ms
t = t_0:dt:T_max;
%% create a stair of 5ms long rectangular current pulses with a gap of
% 10ms and the amplitudes 1uA/cm2, 2uA/cm2, 3uA/cm2, 4uA/cm2, 5uA/cm2
i_stimulus = zeros(1, length(t));
i_stimulus(1:5/dt) = 1; % first stair value 
i_stimulus(5/dt+1:(5/dt+10/dt)) = 0;
i_stimulus((5/dt+10/dt)+1:(5/dt+10/dt)+5/dt) = 2;
i_stimulus((5/dt+10/dt)+5/dt+1:(5/dt+10/dt)+5/dt+10/dt) =0;
i_stimulus((5/dt+10/dt)+5/dt+10/dt+1:(5/dt+10/dt)+5/dt+10/dt+5/dt) = 3;
i_stimulus((5/dt+10/dt)+5/dt+10/dt+5/dt+1:(5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt) = 0;
i_stimulus((5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+1:(5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+5/dt) = 4;
i_stimulus((5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+5/dt+1:(5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+5/dt+10/dt) = 0;
i_stimulus((5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+5/dt+10/dt+1:(5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+5/dt+10/dt+5/dt) = 5;
i_stim1 = i_stimulus;
%% stimulus with 2, 4, 8, 16, 32 microAmps
%% create a stair of 5ms long rectangular current pulses with a gap of
% 10ms and the amplitudes 1uA/cm2, 2uA/cm2, 3uA/cm2, 4uA/cm2, 5uA/cm2
i_stimulus = zeros(1, length(t));
i_stimulus(1:5/dt) = 2; % first stair value 
i_stimulus(5/dt+1:(5/dt+10/dt)) = 0;
i_stimulus((5/dt+10/dt)+1:(5/dt+10/dt)+5/dt) = 4;
i_stimulus((5/dt+10/dt)+5/dt+1:(5/dt+10/dt)+5/dt+10/dt) =0;
i_stimulus((5/dt+10/dt)+5/dt+10/dt+1:(5/dt+10/dt)+5/dt+10/dt+5/dt) = 8;
i_stimulus((5/dt+10/dt)+5/dt+10/dt+5/dt+1:(5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt) = 0;
i_stimulus((5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+1:(5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+5/dt) = 16;
i_stimulus((5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+5/dt+1:(5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+5/dt+10/dt) = 0;
i_stimulus((5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+5/dt+10/dt+1:(5/dt+10/dt)+5/dt+10/dt+5/dt+10/dt+5/dt+10/dt+5/dt) = 64;
i_stim2 = i_stimulus;



%%
% Test model
Temp1 = 6.3; % input temperatures in Celsius
Temp2 = 28;
[V1, t] = HodgkinHuxley(i_stim1, Temp1);
[V2, t] = HodgkinHuxley(i_stim2, Temp2);
%% Plot
% plot model
ref = zeros(size(t));
ref(:) = -70;
head = ['Hodgkin and Huxley Neurons; T_1 = ', num2str(Temp1),' °C and T_2 = ', num2str(Temp2), '°C.'];
plot(t,V1-70, 'r','LineWidth', 1), title(head), ylabel('V(t) in \itmV'), xlabel('t in \itms');
hold on, plot(t,V2-70, 'b','LineWidth', 1)
hold on, plot(t, ref, '--k');
grid on
legend('Action potential (i_{stim}(1))', 'Action potential (i_{stim}(2))', 'V_{rest} = -70 mV')
hold on 

%% plot gating constants 1
m1 = m;
n1 = n;
h1 = h;
figure, title('Gating variables m, n and h for Exp. 1')
hold on
plot(t, m,'LineWidth', 1), ylabel('m, n, h'), xlabel('t in \itms');
hold on, plot(t, n, 'LineWidth', 1)
hold on, plot(t, h, 'LineWidth', 1)
grid on
legend('m', 'n', 'h')

%% plot gating constants 2
m2 = m;
n2 = n;
h2 = h;
figure, title('Gating variables m, n and h for Exp. 2')
hold on
plot(t, m2,'LineWidth', 1), ylabel('m, n, h'), xlabel('t in \itms');
hold on, plot(t, n2, 'LineWidth', 1)
hold on, plot(t, h2, 'LineWidth', 1)
grid on
legend('m', 'n', 'h')

%% plot i_na, i_k over time 2
i_Na2 = i_Na;
i_K2 = i_K;
i_L2 = i_L;
figure, title('HH currents i_{Na} and i_{K} for Exp. 2')
hold on
ylabel('Current I in uA'), xlabel('t in \itms');
hold on, plot(t, i_Na2, 'LineWidth', 1)
hold on, plot(t, i_K2, 'LineWidth', 1)
grid on
legend('i_{Na}', 'i_K')

%% plot i_na, i_k over time 1
i_Na1 = i_Na;
i_K1 = i_K;
i_L1 = i_L;
figure, title('HH currents i_{Na} and i_{K} for Exp. 1')
hold on
ylabel('Current I in uA'), xlabel('t in \itms');
hold on, plot(t, i_Na1, 'LineWidth', 1)
hold on, plot(t, i_K1, 'LineWidth', 1)
grid on
legend('i_{Na}', 'i_K')

%% phase plot 2
figure, title('Phase plot: i_{Na}, i_{K} and i_L against V for Exp. 2')
hold on
ylabel('Current I in \ituA'), xlabel('V in \itmV');
hold on, plot(V, i_Na2, 'LineWidth', 1)
hold on, plot(V, i_K2, 'LineWidth', 1)
hold on, plot(V, i_L2, 'LineWidth', 1)
grid on
legend('i_{Na}', 'i_K', 'i_L')
%xlim([-4 6])

%% phase plot 1
figure, title('Phase plot: i_{Na}, i_{K} and i_L against V for Exp. 1')
hold on
ylabel('Current I in \ituA'), xlabel('V in \itmV');
hold on, plot(V, i_Na1, 'LineWidth', 1)
hold on, plot(V, i_K1, 'LineWidth', 1)
hold on, plot(V, i_L1, 'LineWidth', 1)
grid on
legend('i_{Na}', 'i_K', 'i_L')
%xlim([-4 6])












