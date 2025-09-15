% Alternating Optimization Algorithm Framework:
% Power and SSB periodicity optimization to minimize UE random access delay in LEO satellite system

clear; clc;

% Parameter settings (example values, adjust as needed)
M = 4;              % Number of beams
K = 8;              % Number of cells
U = 50;             % Number of User Equipments (UE)
S = 5;              % Number of epochs
N_iter = 3;        % Maximum iterations
N_MC = 100;         % Number of Monte Carlo samples
P_total = 10;      % Maximum total power of satellite (arbitrary units)
N = 10;             % Maximum SSB periodicity (integer upper bound)

% Initialize power allocation and SSB periodicity randomly or evenly
P = repmat(P_total/M, M, S);     % Initial power equally divided per beam and epoch
theta = randi([1 N], K, S);      % Random initial SSB periodicity for each cell and epoch

% Main optimization loop: alternating between optimizing theta and P
for iter = 1:N_iter
    fprintf('Iteration %d:\n', iter);
    
    % Step 1: Fix power P, optimize SSB periodicity theta
    fprintf('optimize_theta %d:\n', iter);
    theta = optimize_theta(P, theta, M, K, U, S, N_MC, N);
    
    % Step 2: Fix SSB periodicity theta, optimize power P
    fprintf('optimize_power %d:\n', iter);
    P = optimize_power(P, theta, M, K, U, S, N_MC, P_total);
    
    % Optional: Add convergence check and break if converged
end

% Display the resulting optimal power and periodicity after iterations
disp('Final Power Allocation (P):');
disp(P);
disp('Final SSB Periodicity (theta):');
disp(theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to optimize SSB periodicity theta when power is fixed
function theta_opt = optimize_theta(P, theta, M, K, U, S, N_MC, N)
    theta_opt = theta;
    % Loop over each epoch and each cell to find the best periodicity
    for s = 1:S
        for k = 1:K
            best_delay = inf;
            best_theta = theta(k,s);
            % Try all possible periodicity values from 1 to N
            for th = 1:N
                theta_opt(k,s) = th;
                % Estimate expected delay with current theta and fixed P
                fprintf('estimate delay s: %d, k: %d, n: %d\n', s, k, th);
                est_delay = monte_carlo_delay(P, theta_opt, M, K, U, S, N_MC);
                if est_delay < best_delay
                    best_delay = est_delay;
                    best_theta = th;
                end
            end
            theta_opt(k,s) = best_theta; % Update with best periodicity found
        end
    end
end

% Function to optimize power allocation P when SSB periodicity theta is fixed
function P_opt = optimize_power(P, theta, M, K, U, S, N_MC, P_total)
    % P: current power allocation matrix (M x S)
    % theta: current SSB periodicity matrix (K x S)
    % M: number of beams
    % K: number of cells
    % S: number of epochs
    % P_total: total power constraint
    %
    % Returns:
    % P_opt: optimized power allocation (M x S) obeying power constraint
    
    P_opt = zeros(M, S);  % Initialize output

    for s = 1:S
        % Example heuristic: allocate more power to beams serving "difficult" cells
        % Here, 'difficult' is simulated by summing the SSB period for each beam
        % You should map cells to beams as necessary for your scenario
        % For simplicity, assume each beam handles K/M cells (can be customized)
        % Compute an importance score for each beam in this epoch:
        cell_idx = reshape(1:K, M, []);
        beam_scores = zeros(M, 1);
        for m = 1:M
            cells_for_beam = cell_idx(m, :);
            beam_scores(m) = sum(theta(cells_for_beam, s)); % total SSB period in assigned cells
        end
        % Prevent zero scores (avoid divide by zero in normalization)
        beam_scores = beam_scores + 1e-6;

        % Allocate power proportional to beam_scores, keeping total within P_total
        P_opt(:, s) = P_total * (beam_scores / sum(beam_scores));

        % Apply minimum power as needed (e.g., for feasibility)
        min_power = 0; % set a minimum if needed
        for m = 1:M
            if P_opt(m, s) < min_power
                P_opt(m, s) = min_power;
            end
        end

        % Re-normalize if enforced minimum changes total sum
        if sum(P_opt(:, s)) > P_total
            P_opt(:, s) = P_total * (P_opt(:, s) / sum(P_opt(:, s)));
        end
    end
end


% Monte Carlo simulation function to estimate the average UE cell searching delay
% given power allocation P and SSB periodicity theta
function delay_est = monte_carlo_delay(P, theta, M, K, U, S, N_MC)
    % Parameters (examples; replace with your scenario)
    P_th = 0.5;             % SSB detection threshold (customize)
    lambda = 0.03;          % Wavelength [m] example
    Gmax = 10;              % Max antenna gain [linear]
    theta3dB = pi/6;        % Half-power beamwidth [rad]
    % Shadowed-Rician params (choose suitable values)
    m = 2; b = 1; Omega = 3; 

    % Assume: users are randomly assigned to a cell/beams per simulation
    delay_sum = 0;
    for mc = 1:N_MC
        mc_delay = 0;
        for u = 1:U
            % Randomly select cell k and beam m for each UE
            k = randi(K); %%
            m = randi(M);
            s = randi(S);

            % Calculate SSB period for user's cell
            theta_k = theta(k,s);

            % Free space path loss for cell k
            dk = 500e3 + rand*1500e3; % Distance random in LEO [m] %%
            Lk = (lambda/(4*pi*dk))^2;

            % Shadowed-Rician fading sample for cell k
            hk = shadowed_rician(m, b, Omega); % <-- See function below

            % Random UE location relative to beam center for antenna pattern
            theta_m_u = rand * pi; % [0, pi] uniform
            mu = 2.07123 * sin(theta_m_u)/sin(theta3dB); %%
            G_theta = Gmax * ((besselj(1,mu)/(2*mu) + 36*besselj(3,mu)/(mu^3))^2);

            % Calculate received power at selected slot
            P_m_s = P(m,s);
            P_rcv = P_m_s * Lk * hk * G_theta;

            % (3.7) Initial waiting time: uniform [0, theta_k]
            beta_u = rand * theta_k;

            % (3.10) Failure probability for SSB reception
            failure_prob = shadowed_rician_cdf(P_th / (P_m_s * Lk * G_theta), m, b, Omega);

            % (3.9) Number of SSB failures before success (geometric distribution)
            Qu = geornd(failure_prob);

            % (3.8) Additional delay
            gamma_u = Qu * theta_k;

            % (3.6) Cell searching delay
            alpha_u = beta_u + gamma_u;

            mc_delay = mc_delay + alpha_u;
        end
        delay_sum = delay_sum + mc_delay/U; % Average per UE
    end
    delay_est = delay_sum / N_MC;
end

% ----------------- Shadowed-Rician fading random variable generator -----------------
function hk = shadowed_rician(m, b, Omega)
    % This is a stub. Use a custom generator or approximate with Nakagami, Rician, Rayleigh if necessary.
    % For prototyping, use a positive random variable (e.g., gamma) as placeholder:
    hk = gamrnd(m, Omega/(m*b)); 
end

% ----------------- Shadowed-Rician CDF stub (equation (3.2)) -----------------
function cdf_val = shadowed_rician_cdf(x, m, b, Omega)
    % Shadowed-Rician CDF numeric evaluation based on eqn (3.2)
    % Inputs:
    %   x: value where CDF is evaluated (scalar or vector)
    %   m, b, Omega: channel fading parameters
    %
    % Output:
    %   cdf_val: evaluated CDF value(s) at x

    % Parameters
    max_iter = 100;          % Max number of terms in infinite sum
    tol = 1e-8;              % Convergence tolerance

    % Precompute constants
    K = (2*b*m/(2*b*m + Omega))^m / (2*b);
    delta = (Omega / (2*b*m + Omega)) / (2*b);

    % Initialize output
    cdf_val = zeros(size(x));

    % Compute series for each element in x
    for idx = 1:numel(x)
        xi = x(idx);
        sum_val = 0;

        for n = 0:max_iter
            % Compute Pochhammer symbol (m)_n using gamma function
            poch = gamma(m + n) / gamma(m);

            % Coefficient term
            coeff = poch * delta^n * (2*b)^(1+n) / (factorial(n)^2);

            % Lower incomplete gamma function gammainc in MATLAB uses normalized version
            % gammainc(z,s,'lower') is normalized by gamma(s), so multiply by gamma(s)
            lower_gamma = gammainc(xi/(2*b), 1 + n, 'lower') * gamma(1 + n);

            term = coeff * lower_gamma;
            sum_val = sum_val + term;

            % Check convergence
            if abs(term) < tol * abs(sum_val)
                break;
            end
        end

        cdf_val(idx) = K * sum_val;
    end

    % Clip output to [0,1] range for numerical safety
    cdf_val = min(max(cdf_val, 0), 1);
end