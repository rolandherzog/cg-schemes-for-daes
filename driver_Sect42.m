% This is the driver routine for the coupled one-dimensional quasilinear
% heat problem on two subdomains from Section 4.2 of 
% R. Altmann, R. Herzog: "Continuous Galerkin Schemes for Semi-Explicit 
% Differential-Algebraic Equations" (IMA Journal of Numerical Analysis)
clear all

% Set the exponents for the nonlinearities on both subdomains
exponent1 = 3;
exponent2 = 1;

% Setup the 1D heat equation example
f = @(t,x) f_quasilinear_heat_1D(t,x,exponent1,exponent2);
g = @(t,x) g_quasilinear_heat_1D_resistance_and_Dirichlet_onesided(t,x,exponent1,exponent2);

% Setup the initial condition
x0 = zeros(2*41,1);
x0(1:10) = linspace(1,0,10);

% Setup the time span
tspan = linspace(0,0.5,81);

% Set some options
options.degree = 1;

% Solve the problem
solution = daepg(f,g,[],tspan,x0,options);

%% Plot the solution (state) at several time instances
figure(1); clf; hold on
x_plot1 = linspace(0,1,length(x0)/2);
x_plot2 = linspace(1,2,length(x0)/2);
%
for step=8:8:length(tspan)-1
	grid on
    y_plot = solution.eval_state(tspan(step+1));
	y_plot1 = y_plot(1:end/2,:);
	y_plot2 = y_plot(end/2+1:end,:);
 	plot(x_plot1,y_plot1,'LineWidth',2);
	plot(x_plot2,y_plot2,'LineWidth',2);
	axis([0 2 -0.2 1.2]);
	title(sprintf('Solution at time interval %2d',step));
	xlabel('spatial position x');
	drawnow
end

% Evaluate the constraint and multiplier at suitable times
t_plot = solution.times_multiplier;
g_plot = g(0,x0);
g_plot = zeros(size(g_plot,1),length(t_plot));
for i=1:length(t_plot)
	g_plot(:,i) = g(t_plot(i),solution.eval_state(t_plot(i)));
end
lambda_plot = reshape(solution.coeffs_multiplier,size(solution.coeffs_multiplier,1),[]);

%% Plot the associated multiplier - Dirichlet constraint on the left
figure(3); clf; hold on;
plot(t_plot,lambda_plot(1,:),'k.');
grid on
title('Dirichlet constraint multiplier over time');
xlabel('time t');

%% Plot the associated multiplier - first resistance constraint 
figure(5); clf; hold on;
plot(t_plot,lambda_plot(2,:),'k.');
grid on
title('First resistance constraint multiplier over time');
xlabel('time t');

%% Plot the associated multiplier - second resistance constraint
figure(6); clf; hold on;
plot(t_plot,lambda_plot(3,:),'k.');
grid on
title('Second resistance constraint multiplier over time');
xlabel('time t');
