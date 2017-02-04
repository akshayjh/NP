function potential_field
%% This function calculates and plots the potential field for a 50 microm by 50 microm
%  slice in a distance of 10 microm from the point-source

%% def. parameters
rho_medium = 300; % in Ohm.cm
I = 1; % current of point-source in mA
d = 10; % distance of square slice in micro.m 

%% create centered mesh 
[gridx gridy] = meshgrid(0:24, 0:24);
gridx = [fliplr(gridx) gridx];
gridx = [gridx; gridx];
gridy = [flipud(gridy); gridy];
gridy = [gridy gridy];

%% calculate mesh with distances r
gridr = (((gridx.^2 + gridy.^2).^(1/2)).^2 + d^2).^(1/2);

%% calculate potential field
gridphi = phi(rho_medium, I, gridr);


%% plot
figure
surf(gridphi)
title('Potential field', 'FontSize', 24)
xlabel('x in \mum', 'FontSize', 24)
ylabel('y in \mum', 'FontSize', 24)
zlabel('Value of potential \phi in V', 'FontSize', 24)
grid on
colorbar
end