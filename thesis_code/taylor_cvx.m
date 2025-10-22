% 假設 K, C, P_total, N, M 已定義
K = 40;
C = 0.226;
P_total = 100;
N = 600;
M = 4;
cvx_begin
    variables P(K) theta(K)
    expression obj
    obj = sum((1/2 + C * inv_pos(P - C)) .* theta);
    minimize(obj)
    subject to
        sum(P ./ theta) <= P_total;
        P - C >= 1e-6;        % 保證 P > C
        0 < theta <= N;
        sum(inv_pos(theta)) <= M;
cvx_end

