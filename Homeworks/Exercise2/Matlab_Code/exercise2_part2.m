%% Exercise 2
%% Part 2 - differential eq for simple cell model
clear all, clc;
t = linspace(-5, 5, 30); %time Interval
v = linspace(-5, 5, 30); %voltage Interval
[t_grid, v_grid] = meshgrid (t, v); % create a mesh
iso_values = [-2 -1 0 1 2];

% Function handles:
dgl_rhs = @(V, t, R_l, C, Imax) -1/C * (V/R_l + Imax .* sin(t)); % Define the DGL right hand side 
dgl2_rhs = @(V, t, R_l, C, Imax, D) -1/C * (V/R_l + Imax .* sin(t) + D); % Define the DGL right hand side 

isoclines = @(iso, t, R_l, C, Imax, D) -R_l*(Imax * sin(t) + D + C*iso);


%% R_l = 1, C = 1, Imax = 0;

lng = sqrt(dgl_rhs(v_grid , t_grid, 1, 1, 0).^2 + 1); %Length of the vector 
dt = 1 ./ lng ; %Get the horizontal component
dv = dgl_rhs(v_grid , t_grid, 1, 1, 0) ./ lng ; %Get the vertical component
figure, quiver(t, v, dt, dv, 'r'); % Plot the vectors
hold on
title('Slope field for R_l= 1, C= 1, I_{max}= 0'), xlabel('time t in s'), ylabel('Voltage V in V');
xlim([-5, 5]), ylim([-5, 5])
%hold off

% calclulate and plot isoclines:
for i=1:5
    iso_v1(i,:) = isoclines(iso_values(i), t, 1, 1, 0, 0);
end
hold on, plot(t, iso_v1, 'b');
legend('Slope field', 'Isoclines');

%% R_l = 1, C = 1, Imax = 1;

lng = sqrt(dgl_rhs(v_grid , t_grid, 1, 1, 1).^2 + 1); %Length of the vector 
dt = 1 ./ lng ; %Get the horizontal component
dv = dgl_rhs(v_grid , t_grid, 1, 1, 1) ./ lng ; %Get the vertical component
figure, quiver(t, v, dt, dv, 'r');
hold on
title('Slope field for R_l= 1, C= 1, I_{max}= 1'), xlabel('time t in s'), ylabel('Voltage V in V');
xlim([-5, 5]), ylim([-5 5])
%hold off

% calclulate and plot isoclines:
for i=1:5
    iso_v1(i,:) = isoclines(iso_values(i), t, 1, 1, 1, 0);
end
hold on, plot(t, iso_v1, 'b');
legend('Slope field', 'Isoclines');


%% R_l = 1, C = 1, Imax = 0, D = 2;

lng = sqrt(dgl2_rhs(v_grid , t_grid, 1, 1, 0, 2).^2 + 1); %Length of the vector 
dt = 1 ./ lng ; %Get the horizontal component
dv = dgl2_rhs(v_grid , t_grid, 1, 1, 0, 2) ./ lng ; %Get the vertical component
figure, quiver(t, v, dt, dv, 'r');
hold on
title('Slope field for R_l= 1, C= 1, I_{max}= 0, D= 2'), xlabel('time t in s'), ylabel('Voltage V in V');
xlim([-5, 5]), ylim([-5 5])
%hold off

% calclulate and plot isoclines:
for i=1:5
    iso_v1(i,:) = isoclines(iso_values(i), t, 1, 1, 0, 2);
end
hold on, plot(t, iso_v1, 'b');
legend('Slope field', 'Isoclines');

%% R_l = 1, C = 1, I_ex = 1, D = 2;

lng = sqrt(dgl2_rhs(v_grid , t_grid, 1, 1, 1, 2).^2 + 1); %Length of the vector 
dt = 1 ./ lng ; %Get the horizontal component
dv = dgl2_rhs(v_grid , t_grid, 1, 1, 1, 2) ./ lng ; %Get the vertical component
figure, quiver(t, v, dt, dv, 'r');
hold on
title('Slope field for R_l= 1, C= 1, I_{max}= 1, D= 2'), xlabel('time t in s'), ylabel('Voltage V in V');
xlim([-5, 5]), ylim([-5 5])
%hold off

% calclulate and plot isoclines:
for i=1:5
    iso_v1(i,:) = isoclines(iso_values(i), t, 1, 1, 1, 2);
end
hold on, plot(t, iso_v1, 'b');
legend('Slope field', 'Isoclines');

