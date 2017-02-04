function [V, t] = exponentialEuler(h, N, A, B, V_0, t_0)
%% Exponential Euler
% Input arguments
% h : step length
% N : last time point (i.e. max length)
% DGL : differential equation function handle
% V_0 : Value for V(t_0)
% t_0 : starting point for t
% Output arguments:
% V : Voltage vector V, computed with forward Euler
% t : corresponding time vecotr t

% der Exp. Euler hat bei groﬂen Schrittweiten weiter hinten einen
% konstanten Offset

t = t_0:h:N;
V = zeros(1,length(t)); % pre-allocate voltage vector
V(1) = V_0;

DGL_exp = @(V, t) A(t)*V + B(t); % Define the DGL right hand side based on inputs A and B


%% loop to solve differential equation using the exponential Euler
for i=1:length(t)-1
    V(i+1) = V(i) * exp(A(t(i)) .* h) + (B(t(i))/A(t(i))) * (exp(A(t(i)) .* h) - 1);
end



%% plot result:
%str = ['Exponential Euler with step length h= ', num2str(h),'s.'];
%plot(t,V, 'LineWidth', 2), title(str, 'FontSize', 20), ylabel('V(t)', 'FontSize', 16), xlabel('t', 'FontSize', 16);
%t = linspace(t_0, N, 30); %time Interval
%v = linspace(min(V), max(V), 30); %voltage Interval
%[t_grid, v_grid] = meshgrid (t, v); % create a mesh

% plot the slop field
%lng = sqrt(DGL_exp(v_grid , t_grid).^2 + 1); %Length of the vector
%dt = 1 ./ lng ; %Get the horizontal component
%dv = DGL_exp(v_grid , t_grid) ./ lng ; %Get the vertical component
%hold on, quiver(t, v, dt, dv); % Plot the vectors
%xlim([min(t) max(t)]);
%ylim([min(V) max(V)]);
end
