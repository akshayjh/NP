%% Neuroprothetik - Exercise 3
%% Forward Euler:
DGL = @(V, t) 1 - V - t; % Define the DGL right hand side
N = 1; % length
h = 0.01; % step length
t_0 = -4.5;
t = t_0:h:N;
V = zeros(1,length(t));
V(1) = -4;

% loop to compute the forward Euler
for i=1:length(t)-1
    V(i+1) = V(i) + DGL(V(i),t(i)) .* h;
end

% plot result:
plot(t,V, 'LineWidth', 2), title('Forward Euler'), ylabel('V(t)'), xlabel('t');
t = linspace(t_0, N, 30); %time Interval
v = linspace(min(V), max(V), 30); %voltage Interval
[t_grid, v_grid] = meshgrid (t, v); % create a mesh

lng = sqrt(DGL(v_grid , t_grid).^2 + 1); %Length of the vector
dt = 1 ./ lng ; %Get the horizontal component
dv = DGL(v_grid , t_grid) ./ lng ; %Get the vertical component
hold on, quiver(t, v, dt, dv); % Plot the vectors



%% Heun method:
DGL = @(V, t) 1 - V - t; % Define the DGL right hand side
N = 1; % length
h = 0.1; % step length
t_0 = -4.5;
t = t_0:h:N;
V = zeros(1,length(t));
V(1) = -4;

for i=1:length(t)-1
    A = DGL(V(i), t(i));
    V_tilde = V(i) + A * h;
    B = DGL(V_tilde, t(i+1));
    V(i+1) = V(i) + (A+B)/2 .* h;
end

plot(t,V), title('Heun method'), ylabel('V(t)'), xlabel('t');

%% Exponential Euler

% der Exp. Euler hat bei groﬂen Schrittweiten weiter hinten einen
% konstanten Offset
A = @(t) -1;
B = @(t) 1-t;
DGL_exp = @(V, t) A(t)*V + B(t); % Define the DGL right hand side

N = 10; % length
h = 1; % step length
t_0 = -4.5;
t = t_0:h:N;
V = zeros(1,length(t));
V(1) = -4;

for i=1:length(t)-1
    V(i+1) = V(i) * exp(A(t(i)) .* h) + (B(t(i))/A(t(i))) * (exp(A(t(i)) .* h) - 1);
end


figure, plot(t,V, 'LineWidth', 2), title('Exponential Euler'), ylabel('V(t)'), xlabel('t');


% quiver plot
t = linspace(t_0, N, 30); %time Interval
v = linspace(min(V), max(V), 30); %voltage Interval
[t_grid, v_grid] = meshgrid (t, v); % create a mesh

lng = sqrt(DGL(v_grid , t_grid).^2 + 1); %Length of the vector
dt = 1 ./ lng ; %Get the horizontal component
dv = DGL(v_grid , t_grid) ./ lng ; %Get the vertical component
hold on, quiver(t, v, dt, dv); % Plot the vectors




