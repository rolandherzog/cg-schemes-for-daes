function f = f_pendulum(t,x)
% This is an example of a rhs function representing the energy gradient
% in the example of the index-3 pendulum 

gamma = 9.81;

% Determine the number of x-variables (comes from a first-order formulation)
nx = length(x)/2;

% Evaluate f(t,x) = - nabla E(t,x)
f = zeros(2*nx,1);
f(nx,1) = -gamma;,
f(nx+1:end,1) = -x(nx+1:end);
