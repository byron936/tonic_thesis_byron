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

function cdf_val = shadowed_rician_cdf_derivative(x, m, b, Omega)
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
            lower_gamma_derivative = (xi/(2*b))^n * exp(-1 * xi / (2 * b)) / (2 * b);

            term = coeff * lower_gamma_derivative;
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

m1 = 10.1; b1 = 0.126; Omega1 = 0.835; % average
m2 = 0.739; b2 = 0.063; Omega2 = 8.97e-4; % frequent
m3 = 19.4; b3 = 0.158; Omega3 = 1.29; % infrequent

x = linspace(0, 10, 1000);

y1 = shadowed_rician_cdf(x, m1, b1, Omega1);
y2 = shadowed_rician_cdf(x, m2, b2, Omega2);
y3 = shadowed_rician_cdf(x, m3, b3, Omega3);

y_1 = shadowed_rician_cdf_derivative(x, m1, b1, Omega1);
y_2 = shadowed_rician_cdf_derivative(x, m2, b2, Omega2);
y_3 = shadowed_rician_cdf_derivative(x, m3, b3, Omega3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 40;
% average shadow rician
% a = 0.5
D = readmatrix('population_density.txt')';
P_total = 120;
N = 600;
M = 4;
P = ones(K, 1) * P_total / M;
P_th = 1e-14 / 2.2;
S = 5;
lambda = 0.1;

sat_tbl = readmatrix('satellite_positions_3d.csv'); % [Time_Slot, X_m, Y_m, Z_m]
cell_tbl = readmatrix('cell_center_positions_3d.csv'); % [Index, X(m), Y(m), Z(m)]
satellite_positions = sat_tbl(:, 2:4); % [X, Y, Z]
slots_per_epoch = round(size(satellite_positions,1) / S);

Si = 1;
dk = ones(40, 1);
Lk = ones(40, 1);
for k = 1:K
    t0 = (Si - 1) * slots_per_epoch + 1;
    t1 = Si * slots_per_epoch;
    sat_pos = satellite_positions(round((t0+t1)/2),:); % Use center slot of epoch
    cell_center_pos = cell_tbl(k, 2:4);
    dk(k) = norm(cell_center_pos - sat_pos); %unit: meter
    Lk(k) = (lambda/(4*pi*dk(k)))^2;
end

A = P_th ./ P ./ Lk;
a = mean(A);
C = y1(round(a * 100)) - y_1(round(a * 100)) * a;
E = y_1(round(a * 100)) * P_th ./ Lk;

%{
P = linspace(0, 30, 3000);
y = (P .* C + E(1)) ./ (P .* (1 - C) - E(1));
plot(P, y);
%}

function y = f(x, D, theta, C, E)
    y = D .* theta .* (x .* C + E);
end
function y = h(x, C, E)
    y = x .* (1 - C) - E;
end
function psi = Psi(x, beta, u, D, theta, C, E)
    psi = [-1 * f(x, D, theta, C, E) + beta .* h(x, C, E); -1 + u .* h(x, C, E)];
end
function psi_ = Psi_(x, beta, u, D, theta, C, E)
    psi_ = [-1 .* D .* theta .* C + beta .* (1 - C); u .* (1 - C)];
end

function [P_opt, lambda_history] = solve_sum_of_ratios_MN(D, P0, theta, c, E, P_total, xi, epsilon, tol, max_iter)
    K = length(D);
    % α = [β; u]
    beta0 = zeros(K,1);
    u0 = zeros(K,1);
    for k = 1:K
        numerator = P0(k)*c + E(k);
        denominator = P0(k)*(1-c) - E(k);
        beta0(k) = numerator / denominator;
        u0(k) = 1 / denominator;
    end
    alpha = [beta0; u0];
    lambda_history = [];

    for iter = 1:max_iter
        % Step 1: 解子問題（convex problem，根據MN內部公式
        % min sum_k u_k* (P_k c + E_k) - v_k*(P_k(1-c) - E_k)
        % subject to same constraints as original

        cvx_begin quiet
            variable P_var(K)
            % 變數對應alpha參數
            obj = 0;
            for k = 1:K
                u_k = alpha(k);
                v_k = alpha(K+k);
                obj = obj + D(k)*theta(k)*(u_k*(P_var(k)*c + E(k)) - v_k*(P_var(k)*(1-c) - E(k)));
            end
            minimize(obj)
            subject to
                P_var >= 0
                sum(P_var./theta) <= P_total
                P_var.*(1-c) - E > 1e-6 
        cvx_end
        
        % Step 2: 判斷是否收斂
        % 需計算ψ(alpha)，即MN論文中的等式殘值
        psi = zeros(2*K,1);
        for k = 1:K
            numerator = P_var(k)*c + E(k);
            denominator = P_var(k)*(1-c) - E(k);
            psi(k) = numerator - alpha(K+k)*denominator;
            psi(K+k) = alpha(k)*denominator - 1;
        end
        % 判斷收斂
        if norm(psi) < tol
            break;
        end
        
        % Step 2-1: Modified Newton direction
        psi_prime = MN_jacobian(P_var, alpha, c, E); % 用合適函數計算jacobian
        % 這裡如需更嚴謹，請根據問題具體推導雅可比矩陣，此處略以單位近似
        
        % Backtracking line search
        [alpha_next, p_k, lambda_k, i_k] = mn_step2(@(x) MN_psi(P_var, x, c, E), ...
                                                    @(x) MN_jacobian(P_var, x, c, E), ...
                                                    alpha, xi, epsilon, 100, tol);
        lambda_history = [lambda_history; lambda_k];
        alpha
        alpha = alpha_next
    end
    P_opt = P_var;
end

function val = MN_psi(P, alpha, c, E)
    K = length(P);
    val = zeros(2*K,1);
    for k = 1:K
        numerator = P(k)*c + E(k);
        denominator = P(k)*(1-c) - E(k);
        val(k) = numerator - alpha(K+k)*denominator;
        val(K+k) = alpha(k)*denominator - 1;
    end
end

function J = MN_jacobian(P, alpha, c, E)
    % P: Kx1, 問題主變數，於牛頓步中視作固定
    % alpha: 2Kx1, MN法參數 [beta; u]
    % c, E: scalar/vector
    K = length(P);
    J = zeros(2*K,2*K);
    for k = 1:K
        h = P(k)*(1-c) - E(k); % h_k
        % 設定block (2k-1,2k-1)和(2k,2k)為h, 其它0
        J(2*k-1,2*k-1) = h; % derivative of psi_1 w.r.t beta_k
        J(2*k,2*k)     = h; % derivative of psi_2 w.r.t u_k
    end
end

function [alpha_next, p_k, lambda_k, i_k] = mn_step2(psi, psi_prime, alpha_k, xi, epsilon, max_search, tol)
    % Evaluate psi and its jacobian/jacobi at current iterate
    psi_val = psi(alpha_k);
    psi_grad = psi_prime(alpha_k);

    % If converged, return immediately
    if norm(psi_val) < tol
        alpha_next = alpha_k;
        lambda_k = 0;
        p_k = zeros(size(alpha_k));
        i_k = 0;
        return;
    end

    % Step: Compute modified Newton direction
    % p_k = -inv(psi_prime(alpha_k)) * psi(alpha_k)
    %       (for large or ill-conditioned systems, use linsolve)
    p_k = -psi_grad \ psi_val;  % Equivalent to inv(psi_grad) * psi_val

    % Step: Backtracking line search to find i_k
    for i = 0:max_search
        lambda_k = xi^i;
        alpha_trial = alpha_k + lambda_k * p_k;
        psi_trial = psi(alpha_trial);

        % Check Armijo-type sufficient decrease condition
        if norm(psi_trial) <= (1 - epsilon * lambda_k) * norm(psi_val)
            fprintf("%d-th iteration inequality sucess", i);
            i_k = i;
            alpha_next = alpha_trial;
            return;
        end
    end

    % If not found, use smallest step
    i_k = max_search;
    lambda_k = xi^max_search;
    alpha_next = alpha_k + lambda_k * p_k;
end

for iter = 1:1
    %%%%%%%%% Optimize theta %%%%%%%%%%%%%%%%%
    cvx_begin
        variables theta(K)
        expression obj
        obj = sum(D .* (1/2 + (P .* C + E) ./ (P .* (1-C) - E)) .* theta) / sum(D);
        minimize(obj)
        subject to
            sum(P .* inv_pos(theta)) <= P_total;
            0 < theta <= N;
            sum(inv_pos(theta)) <= M;
    cvx_end
    theta = ceil(theta);

    %%%%%%%%% Optimize P %%%%%%%%%%%%%%%%%
    %{
    lb = zeros(K,1); % P >= 0
    
    % 目標函數
    objective = @(P) sum(D .* (1/2 + (P .* C + E) ./ (P .* (1-C) - E) .* theta)) / sum(D);
    
    % 非線性約束
    nonlcon = @(P) deal([
        sum(P ./ theta) - P_total;              % sum(P ./ theta) <= P_total
        -((P .* C + E) - 1e-8);                 % (P .* C + E) > 0，每維都要 > 0
        -((P .* (1-C) - E) - 1e-8)              % (P .* (1-C) - E) > 0，每維都要 > 0
    ], []);
    
    P0 = P; % 初始猜值
    
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    P = fmincon(objective, P0, [], [], [], [], lb, [], nonlcon, options);
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xi = 0.3; epsilon = 0.2; tol = 1e-6; max_iter = 100;
    [P_opt, lambda_hist] = solve_sum_of_ratios_MN(D, P, theta, C, E, P_total, xi, epsilon, tol, max_iter);
    P = P_opt;
end
