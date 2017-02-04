%% Exercise 4 - Neuroprothetik
%% Hodgkin and Huxley Neuron Model
%% Part 1 - Time constants and steady state values
function timeConstants(do_plot)
%% do_plot = 1 : plot graphs

V = linspace(-100, 100, 1000); % Voltage vector V in mV
T_1 = 6.3; % Temperature 1 in degree Celsius
T_2 = 28; % Temperature 2 in degree Celsius

% Rate equations in 1/s; input vector V in mV
% alpha
alpha_m = (2.5 -0.1.*V)./(exp(2.5 - 0.1.*V) - 1); 
alpha_n = (0.1 - 0.01.*V)./(exp(1-0.1.*V) - 1);
alpha_h = 0.07.*exp(-V/20);
% beta
beta_m = 4*exp(-V/18);
beta_n = 0.125 * exp(-V/80);
beta_h = 1./(exp(3-0.1*V) +1);

% temperature correction factor k(T in degree Celsius)
k = @(T) 3^(0.1*(T-6.3));

% function handles to tau and x_inf
tau = @(alpha, beta, k) 1./(k .* (alpha + beta));
x_inf = @(alpha, beta) alpha./(alpha + beta);

% plot 
if do_plot == 1
    
 
%% T1
figure;
head1 = ['\tau-Plots with T_1 = ', num2str(T_1), ' °C.'];
plot(V, tau(alpha_m, beta_m, k(T_1)));
hold on
plot(V, tau(alpha_n, beta_n, k(T_1)));
hold on
plot(V, tau(alpha_h, beta_h, k(T_1)));
ylabel('\tau_x in \itms');
xlabel('Voltage V in \itmV');
legend('\tau_m(T_1)', '\tau_n(T_1)', '\tau_h(T_1)');
title(head1)
grid on

%% T2
figure;
head2 = ['\tau-Plots with T_2 = ', num2str(T_2), ' °C.'];
plot(V, tau(alpha_m, beta_m, k(T_2)));
hold on
plot(V, tau(alpha_n, beta_n, k(T_2)));
hold on
plot(V, tau(alpha_h, beta_h, k(T_2)));
ylabel('\tau_x in \itms');
xlabel('Voltage V in \itmV');
legend('\tau_m(T_2)', '\tau_n(T_2)', '\tau_h(T_2)');
title(head2)
grid on

%% xinf
figure;
plot(V, x_inf(alpha_m, beta_m));
hold on
plot(V, x_inf(alpha_n, beta_n));
hold on 
plot(V, x_inf(alpha_h, beta_h));
ylabel('x_{\infty}');
xlabel('Voltage V in \itmV');
title('x_{\infty}-plots')
grid on
legend('m_{\infty}', 'n_{\infty}', 'h_{\infty}')


% m
% figure;
% subplot(2,3,1), plot(V, tau(alpha_m, beta_m, k(T_1)));
% ylabel('Tau_m in s');
% hold on
% plot(V, tau(alpha_m, beta_m, k(T_2)));
% legend('tau_m(T_1)', 'tau_m(T_2)');
% hold off
% subplot(2,3,4), plot(V, x_inf(alpha_m, beta_m));
% ylabel('x_inf');
% 
% n
% subplot(2,3,2), plot(V, tau(alpha_n, beta_n, k(T_1)));
% xlabel('V in mV');
% hold on
% plot(V, tau(alpha_n, beta_n, k(T_2)));
% legend('tau_n(T_1)', 'tau_n(T_2)');
% hold off
% subplot(2,3,5), plot(V, x_inf(alpha_n, beta_n));
% xlabel('Voltage V in mV');
% 
% h
% subplot(2,3,3), plot(V, tau(alpha_h, beta_h, k(T_1)));
% hold on
% plot(V, tau(alpha_h, beta_h, k(T_2)));
% legend('tau_h(T_1)', 'tau_h(T_2)');
% hold off
% subplot(2,3,6), plot(V, x_inf(alpha_h, beta_h));
% end
end


