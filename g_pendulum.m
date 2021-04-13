function [g,gx] = g_pendulum(t,x)
% This is the constraint function for the pendulum, i.e., 
% x^2 + y^2 = length^2

len = 1;

% Determine the number of x-variables (comes from a first-order formulation)
nx = length(x)/2;

% Evaluate the constraint value
g = sum(x(1:nx).^2) - len^2;

% Evaluate the constraint Jacobian
if (nargout > 1)
	gx = zeros(1,2*nx);
    gx(1,1:nx) = 2.*x(1:nx);
end

