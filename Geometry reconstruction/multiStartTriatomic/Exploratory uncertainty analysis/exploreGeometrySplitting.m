function exploreGeometrySplitting(g)

figure;
scatter(g(:,1), g(:,2));
title('All');
axis([100 300 100 300]);

figure;
components = {'Op_x', 'Op_y', 'Op_z', 'Cp_x', 'Cp_y', 'Cp_z', 'Sp_x', 'Sp_y', 'Sp_z'};
for N=1:9
    subplot(3,3,N);
    indices = logical(repmat([zeros(1,2^(N-1)) ones(1,2^(N-1))], 1, 2^(9-N)));
    scatter(g(indices,1), g(indices,2));
    title(components{N});
    axis([100 200 100 300]);
    box on;
end
end

