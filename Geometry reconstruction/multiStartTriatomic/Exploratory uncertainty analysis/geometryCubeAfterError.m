r12 = linspace(1.15e-10, 1.25e-10, 10);
r23 = linspace(1.56e-10, 1.66e-10, 10);
theta = linspace(170, 180, 10);
[X,Y,Z] = meshgrid(r12, r23, theta);
g = [X(:) Y(:) Z(:)];
p = simulateMomenta(g, [16 12 32], [2 2 2], true);
pp = p(:,4:12);

error = 0.01;
lowerBound = 1-error;
upperBound = 1+error;
col_num = 1;

pp_low = [pp(:, 1:col_num-1) pp(:, col_num)*lowerBound pp(:, col_num+1:end)];
pp_high = [pp(:, 1:col_num-1) pp(:, col_num)*upperBound pp(:, col_num+1:end)];
pp_error = [pp_low; pp_high];