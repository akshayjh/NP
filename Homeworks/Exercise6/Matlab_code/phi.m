function phi = phi(rho, I, r)
%% This function calculates the potential phi at a distance r from a current point-source
%% output : value of potential phi
%% inputs  
%  rho : specific resistance of the medium
%  I   : current of the point source
%  r   : radial distance r

phi = rho/(4*pi) * I ./ r;

end