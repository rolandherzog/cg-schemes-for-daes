% This is the driver routine for the "Simple Circuit Example" from Section 4.1 of
% R. Altmann, R. Herzog: "Continuous Galerkin Schemes for Semi-Explicit 
% Differential-Algebraic Equations" (IMA Journal of Numerical Analysis)
clear all

%% Setup the example
f = @(t,x) f_circuit(t,x);
g = @(t,x) g_circuit(t,x);
x0 = [0; 0];
T = 1.0;

% Exact solution from PhD thesis of Bächle
q2 = @(s) (100*cos(100*s) + 20000*sin(100*s) - 100*exp(-s/2) )/40001;
ref_states = @(s) [sin(100.*s)-q2(s); q2(s)];
ref_final_states = ref_states(T); 
ref_final_lambda = (-2000100*cos(100*T) - 50001*sin(100*T) + 50*exp(-T/2) )/40001;

% Identify the number of constraints
g0 = g(0,x0);
m = length(g0);

% Setup the time step sizes
steps = 2.^[4:10];


%% Loop over polynomial degree
maxDegree = 5;
for deg = 1:maxDegree
    options.degree = deg;
    options.points = 'uniform';
    %options.points = 'Gauss-Lobatto';

    % Prepare a matrix for the final state value for each grid
    % and the total multiplier action
    final_states = zeros(length(x0),length(steps));
    final_multipliers = zeros(m,length(steps));
    max_error_states = zeros(1,length(steps));
    Deltas = zeros(1,length(steps));
    
    % Repeatedly solve the problem on 
    for nstep = 1:length(steps)

        % Determine the time grid for this round
        tspan = linspace(0,T,steps(nstep)+1);
        Deltas(nstep) = tspan(2) - tspan(1);
        
        % Solve the problem on the current time grid
        solution = daepg(f,g,[],tspan,x0,options);

        % Remember the final value of the state vector
        final_states(:,nstep) = solution.eval_state(tspan(end));

        % Remember the final value of the multiplier vector
        % sum over r (corresponds to <lambda,1_I>)
        lam_temp = sum(solution.coeffs_multiplier,2);
        final_multipliers(:,nstep) = lam_temp(:,:,end);

        % reference: integral of lamda over [T-Delta, T]
        T2 = T - Deltas(nstep);
        ref_subIntegral_lambda(nstep) = (-20001*sin(100*T) + 500.01*cos(100*T) - 100*exp(-T/2) )/40001 ... 
                                    - (-20001*sin(100*T2) + 500.01*cos(100*T2) - 100*exp(-T2/2) )/40001;
    end % for nstep = 1:length(steps)
        
    % Evaluate the errors in the final state 
    % and the l2-norm over all state components
    errors_states = (final_states - ref_final_states);
    error_norms_states(deg,:) = vecnorm(errors_states);
    
    % print orders
    disp('orders cG (error of state in t=T)');
    b = log(error_norms_states(deg,:));
    a = log(Deltas);
    approxRate = (b(2:end)-b(1:end-1)) ./ (a(2:end)-a(1:end-1))

    % Evaluate the errors in the final multiplier 
    % and the l2-norm over all multiplier components
    errors_multipliers = (final_multipliers - ref_subIntegral_lambda);
    error_norms_multipliers(deg,:) = abs(errors_multipliers);
    
end % for deg = 1:maxDegree


%% Graphically evaluate the convergence order for the state
figure(1); clf, hold on
grid on
xlabel('grid size \Delta');
title('Error norm in the state vector at final time');
set(gca,'XScale','log');
set(gca,'YScale','log');
for deg = 1:maxDegree
    loglog(Deltas,error_norms_states(deg,:),'LineWidth',1.5);
end
legend('r1','r2','r3','r4','r5');


%% Graphically evaluate the convergence order for the multiplier
figure(2); clf; hold on;
grid on;
xlabel('grid size \Delta');
title('Error norm in the total multiplier action');
set(gca,'XScale','log');
set(gca,'YScale','log');
for deg = 1:maxDegree
    loglog(Deltas,error_norms_multipliers(deg,:),'LineWidth',1.5);
end
legend('r1','r2','r3','r4','r5');
