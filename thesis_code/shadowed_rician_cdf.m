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
y1 = shadowed_rician_cdf_derivative(x, m1, b1, Omega1);
y2 = shadowed_rician_cdf_derivative(x, m2, b2, Omega2);
y3 = shadowed_rician_cdf_derivative(x, m3, b3, Omega3);
hold on
plot(x, shadowed_rician_cdf_derivative(x, m1, b1, Omega1));
plot(x, shadowed_rician_cdf_derivative(x, m2, b2, Omega2));
plot(x, shadowed_rician_cdf_derivative(x, m3, b3, Omega3));
legend('average', 'frequent', 'infrequent');