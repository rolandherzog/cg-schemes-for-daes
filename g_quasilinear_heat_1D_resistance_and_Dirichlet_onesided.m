function [g,gx] = g_quasilinear_heat_1D_resistance_and_Dirichlet_onesided(t,x,exponent1,exponent2)
% This is an example of a constraint function representing Dirichlet boundary
% conditions in the 1D heat equation on the left of the interval, as well as
% a heat resistance conditions at the center vertex.


% Determine the number of spatial grid points on [0,1]
nx = length(x) / 2;

% Determine the spatial mesh size
h = 1/(nx-1);

% Set the heat transfer coefficient
alpha = 1e+1;

% Evaluate the constraint value
g1 = x(1) - 1;
g2 = -x(end/2-1).^exponent1 / h + x(end/2+0).^exponent1 / h + alpha * (x(end/2+0) - x(end/2+1));
g3 =  x(end/2+1).^exponent2 / h - x(end/2+2).^exponent2 / h + alpha * (x(end/2+1) - x(end/2+0));
g = [g1; g2; g3];

% Evaluate the constraint Jacobian
if (nargout > 1)
	gx = zeros(length(g),2*nx);
	gx(1,1) = 1;

	% Evaluate the heat conduction coefficients on the interior boundaries
	kappa11 = exponent1 * x(end/2-1).^(exponent1-1);
	kappa12 = exponent1 * x(end/2+0).^(exponent1-1);
	kappa21 = exponent2 * x(end/2+1).^(exponent2-1);
	kappa22 = exponent2 * x(end/2+2).^(exponent2-1);

	% Setup an auxiliary matrix
	aux = [...
		-kappa11/h, kappa12/h + alpha,           - alpha, 0; ...
		         0,           - alpha, kappa21/h + alpha, -kappa22/h; ...
		];

	gx(2:3,end/2-1:end/2+2) = aux;
end

