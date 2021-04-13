function solution = daepg(f,g,J,tspan,x0,options)
%DAEPG Solves a first-order semi-explicit system of differential algebraic 
%   equations (DAE) of the form
%   
%       J x' = f(t,x) + g_x(t,x)^T \lambda
%         0  = g(t,x)
%
%   Inputs:
%          f  a function handle for the forcing function
%          g  a function handle for the constraint function
%          J  a matrix (default: identity)
%      tspan  points in time defining the time grid
%         x0  initial condition at tspan(1)
%    options  (optional) struct of options
%
%   The forcing function f, when called as f(t,x) with scalar t and column vector x, 
%   must return a column vector.
%
%   The constraint function g, when called as g(t,x) with scalar t and column vector x,
%   must return a column vector,x), as well as the Jacobian g_x(t,x).
%
%   The matrix J is typically non-singular and defaults to the identity matrix.
%   For particular index-3 systems, J = [0, I; -I, 0] with an identity matrix I.
%   
%   tspan is a vector of points in time in increasing order.
%
%   x0 is a vector representing the initial conditions. It need not necessarily be 
%   consistent.
%
%   options.degree is the polynomial degree of the state space. 
%   Its default value is 1.
%
%   options.points determines the Lagrange points and is either 'uniform' (default) or
%   'Gauss-Lobatto'.
%
%   Outputs:
%            solution  a solution structure with the following fields, described from
%                      high-level (user-friendly) to low-level. In the following,
%                      nt = length(tspan) and n = length(x0).
%
%          eval_state  eval_state(t) returns the n-by-m value of the state at time(s) t,
%                      where t can be a vector and m = length(t)
%
%        coeffs_state  the n-by-(degree+1)-by-(nt-1) coefficients for the basis of the state
%                      representing the state part of the solution
%
%   coeffs_multiplier  the l-by-degree-by-(nt-1) coefficients for the basis of the multiplier
%                      representing the multiplier part of the solution
%
%    times_multiplier  the degree * (nt-1) points in time where the multiplier Dirac measures
%                      are supported
%
%                 phi  phi(t) returns an m-by-(degree+1) array of evaluations of the basis 
%                      functions for the state at times(s) t, where t can be a column vector 
%                      and m = length(t)
%
%                 psi  psi(t) returns an m-by-(degree) array of evaluations of the basis
%                      functions for the state test space at time(s) t, where t can be a
%                      column vector and m = length(t)
%
%                  xi  xi(f), where f is a scalar- of vector-valued function, evaluates the
%                      application of the Lagrange multiplier basis on [0,1] to f and returns
%                      an nf-by-degree array, where nf is the dimension of the output of f.


% Check and assign the polynomial degree
if (nargin < 6) || (~isfield(options,'degree')) || (isempty(options.degree))
	r = 1;
else
	r = options.degree;
end

% Check and assign the distribution of Lagrange points
if (nargin < 6) || (~isfield(options,'points')) || (isempty(options.points))
	points = 'uniform';
else
	points = options.points;
end

% Check and assign the matrix J
if (isempty(J))
	J = speye(length(x0));
end

% Setup the Lagrange points on the interval [0,1]
if strcmpi(points,'uniform')
	points = linspace(0,1,r+1)';
elseif strcmpi(points,'Gauss-Lobatto')
	points = gausslobatto(r+1);
end

% Determine the system dimensions
n = length(x0);
m = length(g(tspan(1),x0));

% Determine the indices into the vector of unknown coefficients
% on a single interval
ixx = 1:n*r;
ixlambda = n*r+1:n*r+m*r;

% Determine the time step lengths and the number of intervals
dt = diff(tspan);
N = length(dt);

% Prepare the output
solution.coeffs_state = zeros(n,r+1,N);
solution.coeffs_multiplier = zeros(m,r,N);

% Setup the polynomial bases and the relevant matrices
% over the unit interval [0,1]
[D,M,phi,phiprime,psi,xi,xbar_phi,xbar_psi,xbar_xi] = setup_bases(points);

% Store the basis as part of the solution
solution.phi = phi;
solution.psi = psi;
solution.xi = xi;

% Evaluate the basis psi on the Lagrange points of the basis phi
psi_values = psi(xbar_phi)';

% Set the first coefficient on the very first interval.
% Notice that we are assuming a Lagrangian basis for the solution space
% whose first basis element \phi_0 is supported in the left-most point of the
% interval.
xInitial = x0;

% Loop over all time intervals
for step = 1:N

	% Set the start and end point 
	tstart = tspan(step);
	tend = tspan(step+1);

	% Set up the residual function of the nonlinear system for the current time step
	res = @(z) residual(z(ixx),z(ixlambda),n,m,f,g,J,xInitial,tstart,tend,points,r,D,M,psi_values);

	% Solve the nonlinear system for the combined coefficients z = (x,lambda)
	% on the current interval
	fsolve_options = optimoptions(@fsolve,'Display','off','OptimalityTolerance', 1e-13);
	z0 = [repmat(xInitial,r,1); zeros(m*r,1)];
	z = fsolve(res,z0,fsolve_options);

	% Store the solution's coefficients
	solution.coeffs_state(:,1,step) = xInitial;
	solution.coeffs_state(:,2:r+1,step) = reshape(z(ixx),n,r);
	solution.coeffs_multiplier(:,:,step) = reshape(z(ixlambda),m,r);

	% Remember the coefficients at the right-most point as the initial
	% coefficients for the subsequent interval.
	% Notice that we are assuming a Lagrangian basis for the solution space whose 
	% last basis element \phi_r is supported in the left-most point of the interval.
	xInitial = solution.coeffs_state(:,r+1,step);

end % for step = 1:N


% Assign the evaluation routine for the state x as part of the solution
solution.eval_state = @(t) evaluate_state(t,tspan,phi,solution.coeffs_state);

% Assign feasible evaluation times for the multiplier lambda
solution.times_multiplier = tspan(1:end-1) + dt .* xbar_xi;
solution.times_multiplier = solution.times_multiplier(:);


function res = residual(x,lambda,n,m,f,g,J,xInitial,tstart,tend,points,r,D,M,psi_values)
% This function assembles the residual for the current time step as a function of 
% the unknown coefficients x (w.r.t. the basis \phi_0, ..., \phi_r) 
% and unknown coefficients \lambda (w.r.t. the basis \xi_1, ..., \xi_r).

% Reshape the iterates x and lambda
x = [xInitial, reshape(x,n,r)];
lambda = reshape(lambda,m,r);

% Get the length of the current interval [tstart,tend]
dt = tend - tstart;

% Determine the Lagrange points t_0, ..., t_r on the interval [tstart,tend]
t = tstart + points * dt;

% Evaluate the rhs function f at all Lagrange points t_0, ..., t_r
f_values = zeros(n,r+1);
for j=1:r+1
	f_values(:,j) = f(t(j),x(:,j));
end

% Evaluate the constraint function g and its Jacobian
g_values = zeros(m,r+1);
gx_values = zeros(m,n,r+1);
for j=1:r+1
	[g_values(:,j),gx_values(:,:,j)] = g(t(j),x(:,j));
end

% Evaluate the constraint Jacobian transposed times lambda
gxT_lambda_values = zeros(n,r);
for k=1:r
	gxT_lambda_values(:,k) = gx_values(:,:,k+1)' * lambda(:,k);
end

% Assemble the residuals for the state equation for all test basis function
% \psi_1, ..., \psi_r
res1 = zeros(n,r);
for i=1:r
	res1(:,i) = sum(D(i,:) .* (J*x) - dt * M(i,:) .* f_values,2);
	res1(:,i) = res1(:,i) + sum(psi_values(i,2:end) .* gxT_lambda_values,2);
end

% Assemble the residuals for the constraint equation 
res2 = g_values(:,2:end);

% Assemble the combined residual
res = [res1(:); res2(:)];

