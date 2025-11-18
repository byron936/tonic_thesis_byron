S = 5; % Number of epochs
POP_SIZE = 100; 
MAX_GEN = 100;
GA_GEN = 100;
T_max = 100; % Initial temperature for SA
T_min = 1; % Minimum temperature for SA
cooling_rate = 0.95;
M = 4; % Number of beams
P_total = 120; % Total power budget

pop_den = readmatrix('population_density.txt');
sat_tbl = readmatrix('satellite_positions_3d.csv'); % [Time_Slot, X_m, Y_m, Z_m]
cell_tbl = readmatrix('cell_center_positions_3d.csv'); % [Index, X(m), Y(m), Z(m)]

U = sum(pop_den);
K = length(pop_den);
total_slots = size(sat_tbl, 1);    % Number of time slots

satellite_positions = sat_tbl(:, 2:4); % [X, Y, Z]
slots_per_epoch = round(size(satellite_positions,1) / S);

% Function to calculate expected delay alphauSi based on SSB periodicity thetakSi and power PkSi
function total_expected_delay = calculate_alphauSi(Si, K, thetakSi, PkSi, pop_den, satellite_positions, slots_per_epoch, cell_tbl)
    % size of thetakSi: K * 1
    % size of PkSi: K * 1
    % Implement equations 3.7 to 3.11 here for delay calculation
    % This function returns expected delay for all cells and epochs given thetakSi and PkSi
    P_th = 1e-14 / 2.2; 
    % fr1: 0.410GHz~7.125GHz; fr2: 24.25GHz~71GHz
    lambda = 0.1; % fc = 3GHz
    % Gmax = 40; phi3dB = 0.058; % [6]
    m = 10.1; b = 0.126; Omega = 0.835; % average shadowing [5]

    E_alpha = zeros(K, 1);
    failure_prob = zeros(K, 1);
    dk = zeros(K, 1);
    total_expected_delay = 0;

    for k = 1:K
        % Average over slots in this epoch (pick center slot or average all slots)
        t0 = (Si - 1) * slots_per_epoch + 1;
        t1 = Si * slots_per_epoch;
        sat_pos = satellite_positions(round((t0+t1)/2),:); % Use center slot of epoch
        cell_center_pos = cell_tbl(k, 2:4);
        dk(k) = norm(cell_center_pos - sat_pos); %unit: meter
        Lk = (lambda/(4*pi*dk(k)))^2;

        P = PkSi(k);
        theta = thetakSi(k);

        % SSB reception failure probability (expectation, see eqs)
        % P_th / (P * Lk)
        failure_prob(k) = shadowed_rician_cdf(P_th / (P * Lk), m, b, Omega);

        % Expected initial waiting time
        E_beta = 0.5 * theta;

        % Expected number of failures
        E_Q = failure_prob(k) / (1 - failure_prob(k));

        % Expected additional delay
        E_gamma = E_Q * theta;

        % Expected cell searching delay
        E_alpha(k) = E_beta + E_gamma;
        total_expected_delay = total_expected_delay + E_alpha(k) * pop_den(k);
    end
end

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

% Genetic Algorithm: Fix power P, optimize SSB periodicity thetakSi
function [best_thetakSi] = optimize_SSB_periodicity(Si, PkSi, pre_thetakSi, pop_den, satellite_positions, K, M, POP_SIZE, GA_GEN, slots_per_epoch, cell_tbl)
    % Initialize population: each candidate is a (K x S) matrix of integer periodicities
    population = ones(K, POP_SIZE);
    population(:, 1) = pre_thetakSi;
    for i = 2:POP_SIZE
        while sum(1./population(:, i)) < M - 1 || sum(1./population(:, i)) > M
            population(:, i) = randi([1, K / M * 10], K, 1);
        end
    end

    for gen = 1:GA_GEN
        fitness = zeros(POP_SIZE,1);
        for i = 1:POP_SIZE
            thetakSi = population(:,i);
            % Check constraints: thetakSi integer in [1, slots_per_epochM], vector length constraints from 4.1d
            if sum(1./thetakSi) > M
                fitness(i) = 0;
                continue;
            end
            % Compute expected delay
            totalalphaSi = calculate_alphauSi(Si, K, thetakSi, PkSi, pop_den, satellite_positions, slots_per_epoch, cell_tbl);
            fitness(i) = 1 / totalalphaSi; % Fitness = inverse of average delay
        end
        % Select parents via roulette wheel or tournament selection
        if mod(POP_SIZE / 2, 2) == 0 
            num_parents = POP_SIZE / 2;
        else
            num_parents = POP_SIZE / 2 + 1;
        end
        parents_idx = select_parents(fitness, num_parents); 
        % Crossover & mutation to generate offspring applying constraints
        crossover_rate = 0.9; mutation_rate = 0.1;
        offspring = generate_offspring(population, parents_idx, K, crossover_rate, mutation_rate);
        % Calculate fitness for offspring and select next generation population
        population = select_next_generation(population, K, offspring, fitness, Si, PkSi, pop_den, satellite_positions, slots_per_epoch, cell_tbl);
    end
    % Return best individual
    [~, best_idx] = max(fitness);
    best_thetakSi = population(:, best_idx);
end

% roulette wheel
function parents_idx = select_parents(fitness, num_parents)
    % Ensure fitness is positive
    fitness = fitness - min(fitness);
    fitness = fitness + eps;
    
    % Normalize fitness to get probabilities
    prob = fitness / sum(fitness);
    cum_prob = cumsum(prob);
    
    % Select parents
    parents_idx = zeros(num_parents,1);
    for p = 1:num_parents
        r = rand; % Spin roulette
        idx = find(cum_prob >= r, 1, 'first');
        parents_idx(p) = idx;
    end
end

function offspring = generate_offspring(population, parents_idx, K, crossover_rate, mutation_rate)
    POP_SIZE = length(parents_idx);
    offspring = zeros(K, POP_SIZE);

    for i = 1:2:POP_SIZE
        p1 = population(:, parents_idx(i));
        if i+1 > POP_SIZE
            p2 = population(:, parents_idx(i));
        else
            p2 = population(:, parents_idx(i+1));
        end

        % Perform crossover based on crossover_rate
        if rand <= crossover_rate
            cp = randi([1 K-1]);
            child1 = [p1(1:cp); p2(cp+1:K)];
            child2 = [p2(1:cp); p1(cp+1:K)];
        else
            % No crossover: children are copies of parents
            child1 = p1;
            child2 = p2;
        end

        % Mutation based on mutation_rate
        if rand <= mutation_rate
            mut_idx1 = randi(K);
            child1(mut_idx1) = randi([min(population(:)), max(population(:))]);
        end
        if rand <= mutation_rate
            mut_idx2 = randi(K);
            child2(mut_idx2) = randi([min(population(:)), max(population(:))]);
        end

        offspring(:, i) = child1;
        offspring(:, i+1) = child2;
    end
end

function next_population = select_next_generation(population, K, offspring, fitness, Si, PkSi, pop_den, satellite_positions, slots_per_epoch, cell_tbl)
    % 合併原始 population 及 offspring
    all_individuals = [population, offspring];
    % 需同時計算 offspring fitness，這裡假設 offspring_fitness 已經算好
    all_fitness = [fitness; zeros(size(offspring,2),1)];
    % 你應該在 generate_offspring 後另外算 offspring 適應度
    for i = 1:size(offspring,2)
        % 根據你的計算公式
        individual = offspring(:,i);
        totalalphaSi = calculate_alphauSi(Si, K, individual, PkSi, pop_den, satellite_positions, slots_per_epoch, cell_tbl); 
        all_fitness(size(population,2)+i) = 1 / totalalphaSi;
    end
    % 根據適應度排序（由大到小）
    [~, idx_sorted] = sort(all_fitness,'descend');
    % 取前 POP_SIZE 個
    next_population = all_individuals(:, idx_sorted(1:size(population,2)));
end

% Simulated Annealing: Fix SSB periodicity thetakSi, optimize power allocation PkSi
function [best_P] = optimize_power_allocation(Si, thetakSi, pre_PkSi, pop_den, satellite_positions, slots_per_epoch, K, P_total, T_max, T_min, cooling_rate, cell_tbl)
    % Initialize power allocation Pcurrent with feasible solution
    Pcurrent = pre_PkSi; % Equal power initial guess
    current_delay = calculate_alphauSi(Si, K, thetakSi, Pcurrent, pop_den, satellite_positions, slots_per_epoch, cell_tbl);

    T = T_max;
    best_P = Pcurrent;
    best_delay = current_delay;

    while T > T_min
        % Perturb power allocation to get new candidate, keep sum power <= P_total
        Pnew = perturb_power_allocation(Pcurrent, P_total, K, thetakSi);
        new_delay = calculate_alphauSi(Si, K, thetakSi, Pnew, pop_den, satellite_positions, slots_per_epoch, cell_tbl);

        if new_delay < current_delay
            Pcurrent = Pnew;
            current_delay = new_delay;
            if new_delay < best_delay
                best_P = Pnew;
                best_delay = new_delay;
            end
        else
            % Accept worse solution with probability exp(-(new_delay-current_delay)/T)
            if rand() < exp(-(new_delay - current_delay) / T)
                Pcurrent = Pnew;
                current_delay = new_delay;
            end
        end
        T = T * cooling_rate;
    end
end

function P_new = perturb_power_allocation(P_current, P_total, K, thetakSi)
    % Perturb power allocation vector P_current keeping sum <= P_total
    % Inputs:
    %   P_current: current power allocation (K x 1)
    %   P_total: total power budget (scalar)
    
    P_new = P_current;
    
    % Select two random indices to transfer power between
    idx1 = randi(K);
    idx2 = randi(K);
    while idx2 == idx1
        idx2 = randi(K);
    end
    
    % Choose a small perturbation value delta
    perturbation_rate = 0.1;
    delta = perturbation_rate * P_total * (2*rand - 1); % random between -pert_rate*P_total to +pert_rate*P_total
    
    % Attempt to transfer delta power from idx1 to idx2
    if (P_new(idx1) - delta * thetakSi(idx1) >= 0) && (P_new(idx2) + delta * thetakSi(idx2) >= 0)
        P_new(idx1) = P_new(idx1) - delta * thetakSi(idx1);
        P_new(idx2) = P_new(idx2) + delta * thetakSi(idx2);
    else
        % If invalid, do not change powers (or could try smaller delta)
        % Here no change
    end
    
    % Ensure sum stays approximately P_total by normalization (if needed)
    P_ = P_new ./ thetakSi;
    total_power = sum(P_);
    if abs(total_power - P_total) > 1e-6
        P_ = P_ * (P_total / total_power);
        P_new = P_ .* thetakSi;
    end
end
%%%%%
% Main alternating optimization loop
max_iterations = 100;

Si = 1;
PkSi = ones(K, 1) * (P_total / M);
thetakSi = ones(K, 1) * (K / M);

best_delay = 10000;
best_PkSi = ones(K, 1) * (P_total / K);
best_thetakSi = ones(K, 1) * (K / M);
best_iter = 0;

delay_record = zeros(1, max_iterations);
thetakSi_record = zeros(K, max_iterations);
Pksi_record = zeros(K, max_iterations);

PkSi = [10.298298220186378
10.831882235764285
11.833284237859049
13.302335251705541
15.238787419976633
17.642314040889879
20.512509599664178
10.349172281719978
11.111772694869419
12.339124879087482
14.031021764253477
16.187177897315106
18.807229518212914
21.890734574808850
10.606487373463043
11.133148276039543
12.121557285766263
13.571547647969187
15.482874722835465
17.855216044750268
20.688171369191707
11.895476285950719
12.647837956818403
13.858711593301464
15.527892911606324
17.655100294007092
20.239974851880287
23.282080487493470
13.396577907415006
13.915926957518327
14.890613475678204
16.320473015368602
18.205264334847310
20.544669450659573
23.338293692635162
15.923540114902774
16.665107594763811
17.858608741428700
19.503842187893145
21.600530365370528];
thetakSi = ceil([58
41
32
26
22
18
14
21
19
17
15
14
11
9
15
15
14
12
11
10
8
12
11
11
10
9
7
6
10
10
9
8
8
7
5
8
8
7
7
6]);

%PkSi = ones(K, 1) * (P_total / M);
%thetakSi = ones(K, 1) * (K / M);

cur_delay = calculate_alphauSi(Si, K, thetakSi, PkSi, pop_den, satellite_positions, slots_per_epoch, cell_tbl);
avg_delay = cur_delay / U
%{
for iter = 1:max_iterations
    % 1. Fix power P, optimize SSB periodicity thetakSi with GA
    thetakSi = optimize_SSB_periodicity(Si, PkSi, best_thetakSi, pop_den, satellite_positions, K, M, POP_SIZE, GA_GEN, slots_per_epoch, cell_tbl);
    % 2. Fix periodicity thetakSi, optimize power allocation P with SA
    PkSi = optimize_power_allocation(Si, thetakSi, best_PkSi, pop_den, satellite_positions, slots_per_epoch, K, P_total, T_max, T_min, cooling_rate, cell_tbl);
    % Calculate current delay
    cur_delay = calculate_alphauSi(Si, K, thetakSi, PkSi, pop_den, satellite_positions, slots_per_epoch, cell_tbl);
    avg_delay = cur_delay / U;
    fprintf('Iteration %d, Avg Delay: %f\n', iter, avg_delay);
    delay_record(:, iter) = avg_delay;
    thetakSi_record(:, iter) = thetakSi;
    Pksi_record(:, iter) = PkSi;
    % Check convergence
    if best_delay > avg_delay
        best_delay = avg_delay;
        best_thetakSi = thetakSi;
        best_PkSi = PkSi;
        best_iter = iter;
    end
end

plot(delay_record);
xlabel('Iter');
ylabel('Delay');
title('Total Delay through Iterations');

delay_record = round(delay_record, 4);
Pksi_record = round(Pksi_record);

writematrix(delay_record, 'Delay3.csv');
writematrix(thetakSi_record, 'Periodicity3.csv');
writematrix(Pksi_record, 'Power3.csv');
% final optimal thetakSi and P obtained from above loop
%}