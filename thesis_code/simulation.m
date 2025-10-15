% Function to calculate expected delay alphauSi based on SSB periodicity thetakSi and power PkSi
function [total_expected_delay, E_alpha] = calculate_alphauSi(Si, U, thetakSi, PkSi, ue_positions, satellite_positions, slots_per_epoch, cell_tbl)
    % size of thetakSi: K * 1
    % size of PkSi: K * 1
    % Implement equations 3.7 to 3.11 here for delay calculation
    % This function returns expected delay for all cells and epochs given thetakSi and PkSi
    P_th = 4 * 1e-13; 
    % fr1: 0.410GHz~7.125GHz; fr2: 24.25GHz~71GHz
    lambda = 0.1; % fc = 3GHz
    Gmax = 40; phi3dB = 0.058; % [6]
    m = 10.1; b = 0.126; Omega = 0.835; % average shadowing [5]

    E_alpha = zeros(U, 1);
    failure_prob = zeros(U, 1);
    for u = 1:U
        % Randomly assign this UE to a cell and beam (can use actual assignments if available)
        k = ue_positions(u, 1);

        % Average over slots in this epoch (pick center slot or average all slots)
        t0 = (Si - 1) * slots_per_epoch + 1;
        t1 = Si * slots_per_epoch;
        sat_pos = satellite_positions(round((t0+t1)/2),:); % Use center slot of epoch

        ue_pos = ue_positions(u, 2:4);
        cell_center_pos = cell_tbl(k, 2:4);
        dk = norm(ue_pos - sat_pos); %unit: meter
        Lk = (lambda/(4*pi*dk))^2;

        P = PkSi(k);
        theta = thetakSi(k);

        % Antenna gain 
        vec_ue = ue_pos - sat_pos;

        % Vector from satellite to cell center
        vec_cell = cell_center_pos - sat_pos;
        
        % Compute dot product and norms
        dot_product = dot(vec_ue, vec_cell);
        norm_ue = norm(vec_ue);
        norm_cell = norm(vec_cell);
        
        % Calculate boresight angle in radians
        phi_k_u = acos(dot_product / (norm_ue * norm_cell));
        
        mu = 2.07123 * sin(phi_k_u)/sin(phi3dB);
        G_phi = Gmax * ((besselj(1,mu)/(2*mu) + 36*besselj(3,mu)/(mu^3))^2);

        % SSB reception failure probability (expectation, see eqs)
        failure_prob(u) = shadowed_rician_cdf(P_th / (P * Lk * G_phi), m, b, Omega);

        % Expected initial waiting time
        E_beta = 0.5 * theta;

        % Expected number of failures
        E_Q = failure_prob(u) / (1 - failure_prob(u));

        % Expected additional delay
        E_gamma = E_Q * theta;

        % Expected cell searching delay
        E_alpha(u) = E_beta + E_gamma;
    end
    total_expected_delay = sum(E_alpha);
end