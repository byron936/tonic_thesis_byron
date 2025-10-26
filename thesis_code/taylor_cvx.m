% 假設 K, C, P_total, N, M 已定義
K = 40;
% average shadow rician
% a = 0.5
C = -0.0755;
D = readmatrix('population_density.txt')';
P_total = 120;
N = 600;
M = 4;
P = ones(40, 1) * P_total / M;
P_th = 1e-14 / 2;
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

E = 0.6043 * P_th ./ Lk;

cvx_begin
    variables theta(K)
    expression obj
    obj = sum(D .* (1/2 + (P .* C + E) ./ (P .* (1-C) - E) .* theta)) / sum(D);
    minimize(obj)
    subject to
        sum(P .* inv_pos(theta)) <= P_total;
        0 < theta <= N;
        sum(inv_pos(theta)) <= M;
cvx_end