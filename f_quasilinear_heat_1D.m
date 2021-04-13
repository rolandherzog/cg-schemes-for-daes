function f = f_quasilinear_heat_1D(t,x,exponent1,exponent2)
% This is an example of a rhs function representing the quasilinear 1D heat equation
% discretized by a piecewise linear finite element approach on the interval
% [0,1] and [1,2].

% Define the stiffness matrix as persistent
persistent K1 K2 

% Compute sparse stiffness matrix K only once
if isempty(K1)

	% Determine the number of spatial grid points on [0,1]
	nx = length(x) / 2;

	% Determine the spatial mesh size
	h = 1/(nx-1);

	% Setup the stiffness matrix
	e = ones(nx,1);
	K1 = 1/h^2 * spdiags([-e, 2*e, -e],[-1:1],nx,nx);

	K1(1,1) = K1(1,1) / 2;
	K1(end,end) = K1(end,end) / 2;

	K2 = K1;

end

% Evaluate the rhs value
f = - K1 * x(1:end/2).^exponent1;
f = [f; - K2 * x(end/2+1:end).^exponent2];
