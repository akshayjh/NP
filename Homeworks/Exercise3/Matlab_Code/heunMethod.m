function [V, t] = heunMethod(h, N, DGL, V_0, t_0)
%% Heun method:
% Input arguments
% h : step length
% N : last time point (i.e. max length)
% DGL : differential equation function handle
% V_0 : Value for V(t_0)
% t_0 : starting point for t
% Output arguments:
% V : Voltage vector V, computed with forward Euler
% t : corresponding time vecotr t

t = t_0:h:N;
V = zeros(1,length(t)); % pre-allocate voltage vector
V(1) = V_0;

%% loop to solve differential equation with Heun method
for i=1:length(t)-1
    A = DGL(V(i), t(i));
    V_tilde = V(i) + A * h;
    B = DGL(V_tilde, t(i+1));
    V(i+1) = V(i) + (A+B)/2 .* h;
end

%% plot slope field
t_q = linspace(t_0, N, 30); %time Interval
v_q = linspace(min(V), max(V), 30); %voltage Interval
[t_grid, v_grid] = meshgrid (t_q, v_q); % create a mesh

% plot the slop field
lng = sqrt(DGL(v_grid , t_grid).^2 + 1); %Length of the vector
dt = 1 ./ lng ; %Get the horizontal component
dv = DGL(v_grid , t_grid) ./ lng ; %Get the vertical component
hold on, quiver(t_q, v_q, dt, dv, 'Color', 'k'); % Plot the vectors
xlim([min(t) max(t)]);
ylim([min(V) max(V)]);


%% plot result:
str = ['Heun method with step length h= ', num2str(h),'s.'];
plot(t,V, 'LineWidth', 2), title(str, 'FontSize', 20), ylabel('V(t)', 'FontSize', 16), xlabel('t', 'FontSize', 16);
end

