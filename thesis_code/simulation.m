% Alternating Optimization Algorithm Framework:
% Power and SSB periodicity optimization to minimize UE random access delay in LEO satellite system

clear; clc;

% ----- IMPORT UE AND SATELLITE POSITION DATA -----

% Read UE and satellite position data
ue_tbl = readmatrix('ue_positions_3d.csv'); % [UE_ID, Cell_ID, X_km, Y_km, Z_km]
sat_tbl = readmatrix('satellite_positions_3d.csv'); % [Time_Slot, X_km, Y_km, Z_km]

U = size(ue_tbl, 1);               % Number of UEs
total_slots = size(sat_tbl, 1);    % Number of time slots

% Extract coordinates only
ue_positions = ue_tbl(:, 3:5);     % [X, Y, Z]
satellite_positions = sat_tbl(:, 2:4); % [X, Y, Z]

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
for s = 1:S
    for iter = 1:N_iter
        fprintf('Iteration %d:\n', iter);
        
        % Step 1: Fix power P, optimize SSB periodicity theta
        fprintf('optimize_theta %d:\n', iter);
        theta = optimize_theta(P, theta, M, K, U, s, N_MC, N, ue_positions, satellite_positions);
        
        % Step 2: Fix SSB periodicity theta, optimize power P
        fprintf('optimize_power %d:\n', iter);
        P = optimize_power(P, theta, M, K, U, s, N_MC, P_total, ue_positions, satellite_positions);
        
        % Optional: Add convergence check and break if converged
    end
end

% Display the resulting optimal power and periodicity after iterations
disp('Final Power Allocation (P):');
disp(P);
disp('Final SSB Periodicity (theta):');
disp(theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to optimize SSB periodicity theta when power is fixed
function theta_opt = optimize_theta(P, theta, M, K, U, s, N_MC, N, ue_positions, satellite_positions)
    theta_opt = theta;
    % Loop over each epoch and each cell to find the best periodicity
    for k = 1:K
        best_delay = inf;
        best_theta = theta(k,s);
        %TODO: theta selection algorithm
        % % Try all possible periodicity values from 1 to N
        % for th = 1:N
        %     theta_opt(k,s) = th;
        %     % Estimate expected delay with current theta and fixed P
        %     fprintf('estimate delay s: %d, k: %d, n: %d\n', s, k, th);
        %     est_delay = monte_carlo_delay(P, theta_opt, M, K, U, S, N_MC);
        %     if est_delay < best_delay
        %         best_delay = est_delay;
        %         best_theta = th;
        %     end
        % end
        theta_opt(k,s) = best_theta; % Update with best periodicity found
    end
end

% Function to optimize power allocation P when SSB periodicity theta is fixed
function P_opt = optimize_power(P, theta, M, K, U, s, N_MC, P_total, ue_positions, satellite_positions)
    % P: current power allocation matrix (M x S)
    % theta: current SSB periodicity matrix (K x S)
    % M: number of beams
    % K: number of cells
    % S: number of epochs
    % P_total: total power constraint
    %
    % Returns:
    % P_opt: optimized power allocation (M x S) obeying power constraint
    
    P_opt = zeros(M);  % Initialize output

    % Example heuristic: allocate more power to beams serving "difficult" cells
    % Here, 'difficult' is simulated by summing the SSB period for each beam
    % You should map cells to beams as necessary for your scenario
    % For simplicity, assume each beam handles K/M cells (can be customized)
    % Compute an importance score for each beam in this epoch:
    %TODO: beam-cell allocation
    %TODO: power selection algorithm
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


% Monte Carlo simulation function to estimate the average UE cell searching delay
% given power allocation P and SSB periodicity theta
function delay_est = monte_carlo_delay(P, theta, M, K, U, S, N_MC, ue_positions, satellite_positions)
    % P: (M x S) power allocation matrix
    % theta: (K x S) SSB periodicity matrix
    % U: number of UEs
    % K: number of cells
    % S: number of epochs
    % ue_positions: (U x 3) [X,Y,Z] positions
    % satellite_positions: (T x 3) [X,Y,Z] positions (T = S * slots_per_epoch)
    % Uses: equations (3.6)-(3.10) to compute expected delay analytically

    % Simulation parameters
    P_th = 0.5; lambda = 0.03; Gmax = 10; theta3dB = pi/6;
    m = 2; b = 1; Omega = 3;
    slots_per_epoch = size(satellite_positions,1) / S;

    total_expected_delay = 0;
    for s = 1:S
        epoch_expected_delay = 0;
        for u = 1:U
            % Randomly assign this UE to a cell and beam (can use actual assignments if available)
            k = randi(K);
            m_beam = randi(M);

            % Average over slots in this epoch (pick center slot or average all slots)
            t0 = (s-1)*slots_per_epoch + 1;
            t1 = s*slots_per_epoch;
            sat_pos = satellite_positions(round((t0+t1)/2),:); % Use center slot of epoch

            ue_pos = ue_positions(u,:);
            dk = norm(ue_pos - sat_pos) * 1e3; % [m]
            Lk = (lambda/(4*pi*dk))^2;

            theta_k = theta(k,s);
            P_m_s = P(m_beam,s);

            % Antenna gain (average over angle if needed, here pick mean π/2)
            theta_m_u = pi / 2; % Or sample/average as needed
            mu = 2.07123 * sin(theta_m_u)/sin(theta3dB);
            G_theta = Gmax * ((besselj(1,mu)/(2*mu) + 36*besselj(3,mu)/(mu^3))^2);

            % SSB reception failure probability (expectation, see eqs)
            failure_prob = shadowed_rician_cdf(P_th / (P_m_s * Lk * G_theta), m, b, Omega);

            % Expected initial waiting time
            E_beta = 0.5 * theta_k;

            % Expected number of failures
            E_Q = failure_prob / (1 - failure_prob);

            % Expected additional delay
            E_gamma = E_Q * theta_k;

            % Expected cell searching delay
            E_alpha = E_beta + E_gamma;

            epoch_expected_delay = epoch_expected_delay + E_alpha;
        end
        % Average over UEs
        total_expected_delay = total_expected_delay + epoch_expected_delay / U;
    end

    % Average over epochs
    delay_est = total_expected_delay / S;
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