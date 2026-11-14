% Recreates Figure 5C (Right) from Filo et al. (2026):
% Demonstrates sAIF intrinsic stochastic noise attenuation below the Open-Loop
% Poisson limit (CV = 1/sqrt(E[X])) and shows the Filtered Proportional (fP)
% lower bound using exact Linear Noise Approximation (LNA) equations.

clear; clc; close all;

%% 1. Model Parameters (Matched to Figure 5 in the paper)
gamma = 0.1;     % Output degradation rate
alpha = 2.0;     % Maximal production rate of repression actuation
theta = 1.0;     % Sensing rate constant
kappa = 0.05;    % Repression dissociation constant

% Sweep parameters:
mu_vals  = linspace(1.0, 10.0, 25);       % Sets steady-state output E[X] = mu / theta
eta_vals = logspace(-2, 2, 12);           % Sequestration rate (10^-2 to 10^2)
delta_vals = linspace(0.1, 20.0, 15);     % Degradation rate of Z2 for fP controller

%% 2. Open-Loop Analytical Curve (Equation 13)
% Unimolecular Poisson process: CV[X] = 1 / sqrt(E[X])
E_X_grid = linspace(1.0, 10.0, 200);
CV_open_loop = 1 ./ sqrt(E_X_grid);

%% 3. sAIF Controller Noise via Analytical LNA (Equations 60 & 61)
% Actuation function: h(z2) = alpha / (1 + z2/kappa)
% Inverse at steady state: z2_bar = kappa * (alpha / (gamma * E_X) - 1)
% Derivative gain: sigma_2 = -h'(z2_bar) = (alpha / kappa) / (1 + z2_bar/kappa)^2

E_X_sAIF = [];
CV_sAIF  = [];
eta_color = [];

for i = 1:length(eta_vals)
    eta = eta_vals(i);
    for j = 1:length(mu_vals)
        mu = mu_vals(j);
        E_X = mu / theta; % Robust Perfect Adaptation setpoint
        
        % Ensure steady-state feasibility (alpha > gamma * E_X)
        if alpha > gamma * E_X
            z2_bar = kappa * (alpha / (gamma * E_X) - 1);
            sigma_2 = (alpha / kappa) / (1 + z2_bar / kappa)^2;
            
            % Effective filter cutoff frequency (Equation 61)
            omega_0 = (mu / z2_bar) + eta * z2_bar;
            
            % Analytical CV^2 from LNA (Equation 60)
            num = sigma_2 * theta * (sigma_2 - omega_0) + eta * z2_bar * sigma_2 * (omega_0 + gamma + theta);
            den = (gamma + omega_0) * (gamma * omega_0 + sigma_2 * theta) - eta * z2_bar * sigma_2 * theta;
            
            CV2 = (1 / E_X) * (1 + num / den);
            
            if CV2 > 0
                E_X_sAIF(end+1) = E_X;
                CV_sAIF(end+1)  = sqrt(CV2);
                eta_color(end+1) = log10(eta);
            end
        end
    end
end

%% 4. Filtered Proportional (fP) Controller Noise via LNA (Equation 55)
% Demonstrates the lower noise bound of sAIF when eta -> 0
E_X_fP = [];
CV_fP  = [];

for i = 1:length(delta_vals)
    delta = delta_vals(i);
    % Solve steady-state balance numerically for each delta:
    % alpha / (1 + (theta * E_X / (delta * kappa))) = gamma * E_X
    for E_X_target = linspace(1.0, 10.0, 25)
        % Match effective sensing/degradation ratio to sweep achievable E[X]
        z2_bar = (theta * E_X_target) / delta;
        if z2_bar > 0
            % Check valid parameter regime where h(z2_bar) = gamma * E_X
            alpha_eff = gamma * E_X_target * (1 + z2_bar / kappa);
            if alpha_eff <= alpha * 1.5 % Keep within comparable actuation regime
                sigma_2 = (alpha_eff / kappa) / (1 + z2_bar / kappa)^2;
                
                % Analytical CV^2 for fP controller (Equation 55)
                num_fp = sigma_2 * theta * (sigma_2 - delta);
                den_fp = (gamma + delta) * (gamma * delta + sigma_2 * theta);
                CV2_fp = (1 / E_X_target) * (1 + num_fp / den_fp);
                
                if CV2_fp > 0
                    E_X_fP(end+1) = E_X_target;
                    CV_fP(end+1)  = sqrt(CV2_fp);
                end
            end
        end
    end
end

%% 5. Plot Results (Replicating Figure 5C Right)
figure('Color', 'k', 'Position', [100, 100, 750, 520]);
hold on;

% Plot Filtered P lower bound points
scatter(E_X_fP, CV_fP, 35, [0.65, 0.55, 0.78], 'filled', 'MarkerFaceAlpha', 0.6, ...
    'DisplayName', 'Filtered P (Lower Bound)');

% Plot sAIF points color-coded by log10(eta)
scatter(E_X_sAIF, CV_sAIF, 45, eta_color, 'filled', 'DisplayName', 'sAIF Controller');

% Plot Open Loop curve
plot(E_X_grid, CV_open_loop, 'w-', 'LineWidth', 2.5, 'DisplayName', 'Open Loop (1/\surdE[X])');

% Formatting & Colorbar
colormap(turbo);
cb = colorbar;
cb.Label.String = 'log_{10}(\eta)  [Sequestration Rate]';
cb.Label.FontSize = 11;
clim([-2, 2]);

xlabel('Stationary Mean Output  E[X_1]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Coefficient of Variation  CV[X_1]', 'FontSize', 12, 'FontWeight', 'bold');
title('sAIF Stochastic Noise Attenuation vs. Open Loop (LNA)', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 10);
xlim([1, 10]);
ylim([0, 1]);
grid on;
box on;
hold off;