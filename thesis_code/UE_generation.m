% Parameters
U = 50;          % number of UEs
K = 8;           % number of cells
S = 5;           % number of epochs
slots_per_epoch = 20; % timeslots per epoch (example)
total_slots = S * slots_per_epoch;

% Hexagonal grid parameters (ground plane)
cell_radius = 25; % km, hex radius (distance from center to vertex)

% 3D satellite orbit altitude (km)
orbit_altitude = 700;  % approx LEO altitude

% Calculate hex grid spacing
dx = 3/2 * cell_radius;
dy = sqrt(3) * cell_radius;

% Estimate rows and columns for K cells
cols = ceil(sqrt(K));
rows = ceil(K / cols);

cell_centers = zeros(K, 3); % preallocate 3D: x,y on ground, z=0
idx = 1;
for r = 0:rows-1
    for c = 0:cols-1
        if idx > K, break; end
        x = c * dx;
        if mod(r, 2) == 1
            x = x + dx/2;
        end
        y = r * dy;
        z = 0; % ground level
        cell_centers(idx, :) = [x, y, z];
        idx = idx + 1;
    end
end

% Assign UEs to cells randomly
ue_cells = randi(K, U, 1);

% Generate random UE positions inside cells (circle radius cell_radius)
ue_positions = zeros(U, 3);
for u = 1:U
    center = cell_centers(ue_cells(u), :);
    angle = rand * 2 * pi;
    r_pos = sqrt(rand) * cell_radius;
    x_pos = center(1) + r_pos * cos(angle);
    y_pos = center(2) + r_pos * sin(angle);
    z_pos = 0; % ground level
    ue_positions(u, :) = [x_pos, y_pos, z_pos];
end

% Generate satellite positions in 3D circular orbit (x,y) + fixed altitude z
orbit_radius = 700; % radius of orbit circle in km (distance from Earth center)
angular_velocity = 2 * pi / 90; % one full orbit per 90 slots (example)

satellite_positions = zeros(total_slots, 3);
for t = 1:total_slots
    angle = angular_velocity * t;
    x_sat = orbit_radius * cos(angle);
    y_sat = orbit_radius * sin(angle);
    z_sat = orbit_altitude; % altitude above ground
    satellite_positions(t, :) = [x_sat, y_sat, z_sat];
end

% Write UE positions to CSV (ID, Cell, X,Y,Z)
ue_data = [(1:U)', ue_cells, ue_positions];
ue_filename = 'ue_positions_3d.csv';
ue_header = {'UE_ID', 'Cell_ID', 'X_km', 'Y_km', 'Z_km'};
writecell(ue_header, ue_filename);
writematrix(ue_data, ue_filename, 'WriteMode', 'append');

% Write satellite positions to CSV (Time_Slot, X,Y,Z)
sat_data = [(1:total_slots)', satellite_positions];
sat_filename = 'satellite_positions_3d.csv';
sat_header = {'Time_Slot', 'X_km', 'Y_km', 'Z_km'};
writecell(sat_header, sat_filename);
writematrix(sat_data, sat_filename, 'WriteMode', 'append');

fprintf('3D position files saved: %s and %s\n', ue_filename, sat_filename);
