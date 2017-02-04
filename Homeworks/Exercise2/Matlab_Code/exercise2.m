%% Exercise 2
%% Part 1 - plot slope fields and isocline
clear all, clc;
t = linspace(-5, 5, 30); %time Interval
v = linspace(-5, 5, 30); %voltage Interval
[t_grid, v_grid] = meshgrid (t, v); % create a mesh
iso = [-2 -1 0 1 2];

% Function handles:
one_rhs = @(V, t) 1 - V - t; % Define the DGL right hand side 
two_rhs = @(V, t) sin(t) - V / 1.50; % Define the DGL right hand side 

%% first equation   
lng = sqrt(one_rhs(v_grid , t_grid).^2 + 1); %Length of the vector 
dt = 1 ./ lng ; %Get the horizontal component
dv = one_rhs(v_grid , t_grid) ./ lng ; %Get the vertical component
figure, quiver(t, v, dt, dv, 'k'); % Plot the vectors

% calculate isocline
for i=1:5
    iso_v(i,:) = 1 - t - iso(i);
end

% plot
hold on
plot(t, iso_v), title('Slope field and isoclines for eq. 1');
xlim([-5,5]), xlabel('Time t in s');
ylim([-5,5]), ylabel('Voltage V in V');
legend('Slope field', 'isocline, dV/dt= -2V/s', 'isocline, dV/dt= -1V/s', 'isocline, dV/dt= 0V/s', 'isocline, dV/dt= 1V/s', 'isocline, dV/dt= 2V/s');
hold off



%% second equation
lng = sqrt(two_rhs(v_grid , t_grid).^2 + 1); %Length of the vector 
dt = 1 ./ lng ; %Get the horizontal component
dv = two_rhs(v_grid , t_grid) ./ lng ; %Get the vertical component
figure, quiver(t, v, dt, dv, 'k'); % Plot the vectors

% calculate isocline
for i=1:5
    iso_v(i,:) = 1.5*(sin(t) - iso(i));
end

% plot
hold on
plot(t, iso_v), title('Slope field and isoclines for eq. 2');
xlim([-5,5]), xlabel('Time t in s');
ylim([-5,5]), ylabel('Voltage V in V');
legend('Slope field', 'isocline, dV/dt= -2V/s', 'isocline, dV/dt= -1V/s', 'isocline, dV/dt= 0V/s', 'isocline, dV/dt= 1V/s', 'isocline, dV/dt= 2V/s');
hold off







