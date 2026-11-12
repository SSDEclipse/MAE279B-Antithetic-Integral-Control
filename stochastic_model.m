% --- Script 2: Gillespie SSA Simulation for sAIF Controller ---
clear; clc;

% Biological Parameters
gamma = 0.1;     % Output degradation rate constant
mu    = 5.0;     % Reference product4ion rate
theta = 1.0;     % Sensing rate constant
eta   = 1.0;     % Sequestration rate constant
alpha = 10.0;    % Max actuation production rate
kappa = 1.0;     % Repression threshold constant

% Expected steady-state output
x_setpoint = mu / theta; 

% Reaction Stoichiometry Matrix (3 species: X, Z1, Z2)
% Reactions: 1: Prod X, 2: Deg X, 3: Prod Z1, 4: Prod Z2, 5: Sequestration
S = [ 1, -1,  0,  0,  0 ;
    0,  0,  1,  0, -1 ;
    0,  0,  0,  1, -1 ];

% Initial condition [X, Z1, Z2]
Y = [10; 5; 5];
t = 0;
t_max = 5000;    % Total simulation time
burn_in = 1000;  % Time to discard for transient burn-in

% Preallocate memory for trajectory storage
max_samples = 1e6;
time_vec = zeros(max_samples, 1);
state_history = zeros(max_samples, 3);

idx = 1;
time_vec(idx) = t;
state_history(idx, :) = Y';

% Main Gillespie Loop
while t < t_max && idx < max_samples
    X  = Y(1);
    Z1 = Y(2);
    Z2 = Y(3);

    % Propensity Functions (Negative Actuation via Repression)
    h_Z2 = alpha / (1 + (Z2 / kappa)); % Repression function

    propensities = [
        h_Z2 ;          % R1: 0 -> X
        gamma * X ;     % R2: X -> 0
        mu ;            % R3: 0 -> Z1
        theta * X ;     % R4: X -> X + Z2
        eta * Z1 * Z2   % R5: Z1 + Z2 -> 0
        ];

    a0 = sum(propensities);
    if a0 == 0, break; end

    % Draw time to next reaction (exponential distribution)
    r1 = rand();
    dt = -log(r1) / a0;

    % Determine which reaction occurs
    r2 = rand() * a0;
    cum_prop = cumsum(propensities);
    rxn_idx = find(cum_prop >= r2, 1, 'first');

    % Update State
    Y = Y + S(:, rxn_idx);
    t = t + dt;

    idx = idx + 1;
    time_vec(idx) = t;
    state_history(idx, :) = Y';
end

% Trim unused array space
time_vec = time_vec(1:idx);
state_history = state_history(1:idx, :);

% Filter transient steady-state data (t > burn_in)
steady_idx = time_vec > burn_in;
X_steady = state_history(steady_idx, 1);

% Statistical Evaluations
mean_X = mean(X_steady);
std_X  = std(X_steady);
CV_X   = std_X / mean_X;

fprintf('Desired Setpoint (mu/theta): %.2f\n', x_setpoint);
fprintf('Simulated Mean E[X]: %.3f\n', mean_X);
fprintf('Simulated Standard Deviation: %.3f\n', std_X);
fprintf('Simulated Coefficient of Variation (CV): %.4f\n', CV_X);

% Plot Trajectory
figure;
plot(time_vec, state_history(:, 1), 'b-', 'LineWidth', 1.2); hold on;
yline(x_setpoint, 'r--', 'Setpoint (\mu/\theta)', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Output Species Concentration (X)');
title('sAIF Stochastic Simulation (Gillespie SSA)');
grid on;