function [D,M,phi,phiprime,psi,xi,xbar_phi,xbar_psi,xbar_xi] = setup_bases(points)
% This function sets up the bases 
%   \phi_1, ..., \phi_{r+1}  for the state solution space (and its derivative),
%   \psi_1, ..., \psi_r      for the state test space,
%   \xi_1, ..., \xi_r        for the multiplier space 
% as well as the associated matrices 
%   D_{ij} = \int_I \phi'_j(t) \psi_i(t) dt 
%   M_{ij} = \int_I \phi_j(t)  \psi_i(t) dt 
% with i = 1, ..., r and j = 1, ..., r+1 over the unit interval I = [0,1].
%
% NOTICE that the 'j' indices (for \phi_j) differ from those in the paper by 1. 
%
% The input 'points' contains the support points for the Lagrangian basis 
% in [0,1]. We assume point(1) = 0.

% Get the degree from the number of Lagrangian points
r = length(points) - 1;

% Prepare the output
D = zeros(r,r+1);
M = zeros(r,r+1);

% Determine the interpolation/evaluation points
xbar_phi = points;
xbar_psi = xbar_phi(2:end);
xbar_xi = xbar_psi;

% Loop over all row indices
for i=1:r

	% Setup the function psi_i
	psi_i = @(x) Lagrangian(x,xbar_psi,i);

	% Loop over all column indices
	for j=1:r+1

		% Setup the function phi_j
		phi_j = @(x) Lagrangian(x,xbar_phi,j);

		% Setup the integrand for the matrix entry in M
		f = @(x) psi_i(x) .* phi_j(x);

		% Numerically evaluate the integral for matrix entry in M
		M(i,j) = integral(f,0,1);


		% Setup the function phi'_j
		phiprime_j = @(x) second(@(x) Lagrangian(x,xbar_phi,j),x);

		% Setup the integrand for the matrix entry in D
		f = @(x) psi_i(x) .* phiprime_j(x);

		% Numerically evaluate the integral for matrix entry in D
		D(i,j) = integral(f,0,1);

	end % for j=1:r+1

end % for i=1:r

% Setup an array of all phi functions
phi = @(x) Lagrangian(x,xbar_phi,1);
for j=2:r+1
	phi = @(x) [phi(x), Lagrangian(x,xbar_phi,j)];
end

% Setup an array of all derivatives of phi functions
phiprime = @(x) second(@(x) Lagrangian(x,xbar_phi,1),x);
for j=2:r+1
	phiprime = @(x) [phiprime(x), second(@(x) Lagrangian(x,xbar_phi,j),x)];
end

% Setup an array of all psi functions
psi = @(x) Lagrangian(x,xbar_psi,1);
for i=2:r
	psi = @(x) [psi(x), Lagrangian(x,xbar_psi,i)];
end

% Setup an array of all xi functionals
xi = @(f) feval(f,xbar_xi(1));
for k=2:r
	xi = @(f) [xi(f), feval(f,xbar_xi(k))];
end

% Helper function returning only the 2nd argument of a function
function val = second(f,x)
[~,val] = f(x);

