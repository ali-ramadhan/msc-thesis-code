function exploration2(col_nums, error)
for col_num = col_nums
    exploration2_single(col_num, error);
end
end


function exploration2_single(col_num, error)

% Choose geometry parameters
r12 = linspace(1.15e-10, 1.25e-10, 10);
r23 = linspace(1.56e-10, 1.66e-10, 10);
theta = linspace(170, 180, 10);

% Create list of geomtries (really a 3-column matrix)
[X,Y,Z] = meshgrid(r12, r23, theta);
g = [X(:) Y(:) Z(:)];

% Coulomb explode them all and only keep the momentum results
p = simulateMomenta(g, [16 12 32], [2 2 2], true);
p = p(:,4:12);

lowerBound = 1-error;
upperBound = 1+error;

% Create the lower and upper bounds for the momentum vectors, changing
% only one component.
p_low = [p(:, 1:col_num-1) p(:, col_num)*lowerBound p(:, col_num+1:end)];
p_high = [p(:, 1:col_num-1) p(:, col_num)*upperBound p(:, col_num+1:end)];
p_error = [p_low; p_high];

% Reconstruct each one now.
masses = [16 12.011 32.06];
charges = [2 2 2];
filenamePrefix = strcat('EXP2_COL', num2str(col_num));
multiStartTriatomic(p_error, masses, charges, filenamePrefix, 1, 1);
end