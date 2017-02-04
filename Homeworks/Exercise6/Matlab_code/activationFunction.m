function activationFunction
%% This function 
%  calculates and plots the external potential, the electric field and the
%  activation function along a 50 micro.m piece of aon positioned 10
%  micro.m from a current point-source. 
%  

%% def. parameters
rho_medium = 300; % in Ohm.cm
I = -1; % current of point-source in mA
d = 10; % distance of square slice in micro.m 

h = 0.001; % step size
x = -25:h:25;
r = (d^2 + x.^2).^(1/2);

%% external potential
v_ext = phi(rho_medium, I, r);

% electric field
e_field = -diff(v_ext)/h;

% activation function
activation_fnct = diff(e_field)/h;


%% create plots
%% external potential
figure, plot(x, v_ext, 'b','LineWidth', 2)
title('a) external potential', 'FontSize', 24), xlabel('x in \mum', 'FontSize', 24), ylabel('Potential \phi in V', 'FontSize', 24);
grid on
set(gca,'FontSize',20);
figure, plot(x(1:end-1), e_field, 'r','LineWidth', 2)
title('b) electric field', 'FontSize', 24), xlabel('x in \mum', 'FontSize', 24), ylabel('Electric field E in V/m', 'FontSize', 24);
grid on
set(gca,'FontSize',20);
figure, plot(x(1:end-2), activation_fnct, 'g','LineWidth', 2)
title('c) activation function', 'FontSize', 24), xlabel('x in \mum', 'FontSize', 24), ylabel('Activation function in V/m^2', 'FontSize', 24);
grid on
set(gca,'FontSize',20);


end