function val = evaluate_state(t,tspan,phi,coeffs)
% This function returns the value of the solution at time(s) t, which is
% described by the coefficients coeffs belonging to the basis
% functions \phi_1, ..., \phi_{r+1} on intervals described by tspan.
%
% NOTICE that the 'j' indices (for \phi_j) differ from those in the paper by 1. 

% Get the dimensions
n = size(coeffs,1);
r = size(coeffs,2) - 2;

% Make t a column vector
t = t(:);

% Determine the time step lengths 
dt = diff(tspan);

% Locate the t values within tspan
[ix,jx] = find(tspan(1:end-1) <= t & t <= tspan(2:end));

% Make sure each t is found only once
[ix,reduction] = unique(ix,'stable');
jx = jx(reduction);

% Convert each evaluation point to [0,1] within the respective interval
normalized_points = (t(ix) - tspan(jx)') ./ dt(jx)';

% Evaluate all basis functions 
basis_values = phi(normalized_points);

% Assign the values
val = nan(n,length(t));
for i=1:length(ix) 
	val(:,ix(i)) = coeffs(:,:,jx(i)) * basis_values(i,:)';
end

