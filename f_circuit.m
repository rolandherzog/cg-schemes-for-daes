function f = f_circuit(t,x)
% This is an example of a constraint function pertaining to a linear circuit example

% Evaluate the rhs value pertaining to the ODE part
f = [-sin(100*t); -x(2) - sin(100*t)];

