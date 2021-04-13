function [g,gx] = g_cstr(t,x)
% This is an example of a constraint function pertaining to a linear circuit example

% Evaluate the constraint value
g = x(1,:) + x(2,:) - sin(100*t);

% Evaluate the constraint Jacobian
if (nargout > 1)
	gx = ones(1,2);
end

