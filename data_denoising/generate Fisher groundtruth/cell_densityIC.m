function cell_dense0 = cell_densityIC(x)
% creates an initial condition cell density over 
% given x values
shift = 1000;
a = -1000;
b = 1;
c = 0.6;
cell_dense0 = a/(c*sqrt(2*pi))*exp(-0.5*(x-b).^2/c.^2) + shift;
end % cell_density