% This is the driver routine for the pendulum example (DAE of index 3) from Section 4.3 of 
% R. Altmann, R. Herzog: "Continuous Galerkin Schemes for Semi-Explicit 
% Differential-Algebraic Equations" (IMA Journal of Numerical Analysis)
clear all

% Setup the example
f = @(t,x) f_pendulum(t,x);
g = @(t,x) g_pendulum(t,x);
x0 = [1;0;0;0];
J = [0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, 0; 0, -1, 0, 0];

% Identify the number of constraints
g0 = g(0,x0);
m = length(g0);

% Setup the time step sizes
steps = 2.^[4:11];


%% Loop over polynomial degree
maxDegree = 3;
%
for deg = 1:maxDegree
    options.degree = deg;    
    options.points = 'uniform';

    % Prepare a matrix for the final state value for each grid
    % and the total multiplier action
    final_states = zeros(length(x0),length(steps));
    action_multipliers = zeros(m,length(steps));
    Deltas = zeros(1,length(steps));

    % Repeatedly solve the problem on 
    for nstep = 1:length(steps)

        % Determine the time grid for this round
        tspan = linspace(0,3,steps(nstep)+1);
        Deltas(nstep) = tspan(2) - tspan(1);

        % Solve the problem on the current time grid
        solution = daepg(f,g,J,tspan,x0,options);

        % Remember the final value of the state vector
        final_states(:,nstep) = solution.eval_state(tspan(end));
    
    end % for nstep = 1:length(steps)

    % Evaluate the errors in the final state 
    % and the l2-norm over all state components
    errors_states = (final_states - final_states(:,end));
    errors_states = errors_states(:,1:end-1);
    error_norms_states(deg,:) = vecnorm(errors_states);
    approxRate = (log(error_norms_states(deg,1:end-1)) - log(error_norms_states(deg,2:end)))./log(2)

end


%% Graphically evaluate the convergence order for the state
figure(1); clf; hold on;
grid on;
xlabel('grid size \Delta');
title('Error norm in the state vector at final time');
for deg = 1:maxDegree
    loglog(Deltas(1:end-1),error_norms_states(deg,:),'LineWidth',1.5);
end
set(gca,'XScale','log');
set(gca,'YScale','log');
legend('r1','r2','r3');
