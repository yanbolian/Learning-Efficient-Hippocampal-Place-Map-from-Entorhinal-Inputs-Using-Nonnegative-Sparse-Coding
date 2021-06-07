%% This cdoe uses sparse coding (LCA) to learn place fields given prefixed entorhinal cells
% Author: Yanbo Lian
% Date: 2021-06-07

clear
close all;

%% Set up the environment

% Size of the environment
x_size = 1; % unit: m
y_size = 1; 

% Number of discrete points that represent the environment
Nx = 32;
Ny = 32;
L = Nx * Ny;

%% Generate entorhinal cells (grid cells or weakly spatial cells) in the network

%----------------------------- Grid cells---------------------------%
% Lower and upper spacing of the grid fields
lambda_lower = 0.28; % 28 cm is used in Solstad et al. 2006
lambda_upper = 1; % 73 cm is used in Solstad et al. 2006
spacing_ratio = 1.42; % ratio for the geometric progression of the grid spacings (Wei et al. 2015)

% All possible spacings in the network
num_lambda = floor(log(lambda_upper/lambda_lower) / log(spacing_ratio)) + 1;
lambdas = zeros(num_lambda, 1);
for i = 1 : num_lambda
    lambdas(i) = lambda_lower * spacing_ratio^(i-1);
end

% All possible orientations in the network
% Orientations is no larger than pi/3 because of the triangular pattern
num_orientation = 6;
orientations = pi/3 * rand([num_orientation, 1]); % Uniform distribution of the orientations

% All possible phase offsets in the network
num_phase_x = 5; % possible phases along x-axis
num_phase_y = 5; % possible phases along y-axis
num_phase = num_phase_x * num_phase_y;

num_grid_cell = num_lambda * num_orientation * num_phase;
G = zeros(L, num_grid_cell); % each column is a grid field for a grid cell

i_grid_cell = 0; % Index of grid cells in the network
for i_lambda = 1 : num_lambda
    lca_lambda = lambdas(i_lambda); % grid spacing (unit: m)
    orientations = pi/3 * rand([num_orientation, 1]); % Randomly sample orientations
    
    % Linearly sample the phase offsets for each lambda (spacing)
    [phases_X, phases_Y] = meshgrid( 0 : lca_lambda/num_phase_x : lca_lambda-lca_lambda/num_phase_x, ...
                                    0 : lca_lambda/num_phase_y : lca_lambda-lca_lambda/num_phase_y);
    
    for i_orientation = 1 : num_orientation
        orientation = orientations(i_orientation); % orientation of the grid pattern (unit: radian)
        for i_phase = 1 : num_phase
            i_grid_cell = i_grid_cell + 1;
            phase = [phases_X(i_phase), phases_Y(i_phase)]; % grid phase - (x0,y0) (unit: m)
            figure_on = 0; % Display the grid field or not
            grid = generate_2D_grid_field(x_size, y_size, Nx, Ny, lca_lambda, orientation, phase, figure_on);
            G(:, i_grid_cell) = reshape(grid, Nx*Ny, 1);    
        end
    end
end

%--------------------- Weakly spatial cells------------------------%
num_weakly_cell = 600;
W = zeros(L, num_weakly_cell); % each column is a grid field for a grid cell
for i_weakly_cell = 1 : num_weakly_cell
    figure_on = 0; % Display the field or not
    weakly = generate_2D_weakly_modulated(x_size, y_size, Nx, Ny, figure_on);
    W(:, i_weakly_cell) = reshape(weakly, Nx*Ny, 1);
end

%% Use grid cells or weakly spatial cells as the entorhinal input
E = G; num_entorhinal_cell = num_grid_cell;
% E = W; num_entorhinal_cell = num_weakly_cell;

figure(1); display_matrix(E,3); title('Entorhinal cells'); colormap(jet_modified);colorbar

%% Create place cells and training parameters
num_place_cell = 100; % Number of place cells in the network

n_epoch = 2e4; 3e4; % Number of epoches; more epoches and smaller learning rate for weakly spatial cells
A_eta = 3e-2; 1e-2;% learning rate of A
lca_lambda = 0.3;
lca_tau = 10; % ms
dt = 0.8; % ms
lca_eta = dt/lca_tau;
lca_n_iter = 200;
display_every = 500; % Frequency of generating plot

%% Definitions of symbols
A = normalize_matrix(randn(num_entorhinal_cell, num_place_cell)); % Connections between input and output layer
X = zeros(L, 1); % Input matrix: each column is an image patch in the vector form (sz*sz, 1)
U = randn(num_place_cell, 1); % Membrane potentials of M neurons
S = rand(num_place_cell, 1); % Firing rates (Response) of M neurons

figure(2); display_matrix(normalize_matrix(A),1); title('A'); colormap(gray);colorbar
figure(3); display_matrix(normalize_matrix(E * A),1); title('EA'); colormap(gray);colorbar

%% main loop
residual = 0;
x = 0;
s = 0;

for i_epoch = 1 : n_epoch       

    % Generate a random spatial location
    temp = zeros(Ny, Nx);
    temp(randi([1, Ny]), randi([1, Nx])) = 1;
    X_data(:,i) = reshape(temp, L, 1); % Take uniform input
       
    % Response of grid cells
    X = E' * X_data;
    
    % Compute the response of hippocampal
    [S, U, S_his, U_his] = sparse_coding_by_LCA(...
            X, A, lca_lambda, lca_eta, lca_n_iter);

    R = X - A * S; % Calculate residual error
    
    % update bases
    dA1 = R * S';
    A = A + A_eta * dA1;
    A = max(A, 0); % A is non-negative
    A = normalize_matrix(A, 'L2 norm', 1); % Normalize each column of the connection matrix
    
    % Save R_average, X_average and A_error
    R_average(i_epoch) = sum(sum(R.^2)) / L;
    X_average(i_epoch) = sum(sum(X.^2)) / L;
    residual = residual + R_average(i_epoch);
    x = x + X_average(i_epoch);
    s = s + sum(S(:)~=0) / num_place_cell;
    
    fprintf('Iteration %6d: Percentage of active units: %3.1f%%, MSE: %3.1f%%\n',...
            i_epoch, 100 * sum(S(:) ~= 0) / num_place_cell, 100 * R_average(i_epoch) / X_average(i_epoch)...
                );
       
	% Display connection matrix, A, and responses, S
    if (mod(i_epoch,display_every) == 0)
        figure(1); display_matrix(E,3); title('Grid cells'); colormap(jet_modified);colorbar
        figure(2); display_matrix(normalize_matrix(A),1); title('A'); colormap(gray);colorbar
        figure(3); display_matrix(normalize_matrix(E * A),1); title('EA'); colormap(gray);colorbar

        fprintf('Percentage of active units: %3.1f%%\n', 100 * s / i_epoch);
        fprintf('MSE: %3.1f%%\n', 100 * residual / x);
        
        pause(0.3)
    end
end

%% Recover the place fields
n_input = 1e5;
X_recover = zeros(L, n_input); % each column is a vector with length L that represents Nx*Ny input location

loc_index = randi([1, L], [n_input, 1]); % The index where the animal is at for each input

for i_input = 1 : n_input
    X_recover(loc_index(i_input), i_input) = 1;
end

X = E' * X_recover; % Response of grid cells
[S, U] = sparse_coding_by_LCA(X, A, lca_lambda, lca_thresh_type, lca_eta, lca_n_iter); % Compute the response of place cells
place_field_recover = X_recover(:, 1:n_input) * S' ./ repmat(sum(S,2)', L, 1); % Using reverse correlation / STA to recover the place field of cells

figure(4);
display_matrix(normalize_matrix(place_field_recover), 3);
colormap(jet_modified); colorbar
title('Place field (learned by sparse coding)'); 