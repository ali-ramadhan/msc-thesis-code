function exploration1plot(geometries)
% r_CO=131.677 pm, r_CS=187.161 pm, theta=169.42 deg
r12 = geometries(:,2);
r23 = geometries(:,3);
theta = geometries(:,4);

subplot(3,2,1);
histogram(r12);
xlabel('r_{CO} (pm)');
ylabel('counts');

subplot(3,2,3);
histogram(r23);
xlabel('r_{CS} (pm)');
ylabel('counts');

subplot(3,2,5);
histogram(theta);
xlabel('\theta_{OCS} (deg)');
ylabel('counts');

subplot(3,2,2);
s = scatter(r12, r23, 'o');
s.MarkerEdgeColor = [0.000, 0.447, 0.741];
s.MarkerEdgeAlpha = 0.6;
xlabel('r_{CO} (pm)');
ylabel('r_{CS} (pm)');
axis tight;
box on;

hold on;
[k,v] = convhull(r12, r23);
p = plot(r12(k), r23(k), 'r:', 'LineWidth', 1.5);
p.Color = [0.635, 0.078, 0.184, 0.6];
fprintf('r12*r23 convex hull area: %f [pm^2]\n', v);

hold on;
shp = alphaShape(r12, r23, 15);
bf = boundaryFacets(shp);
p = plot(shp.Points(bf,1), shp.Points(bf,2), '-', 'LineWidth', 1.5);
p.Color = [0.494, 0.184, 0.556, 0.6];
fprintf('r12*r23 alpha shape area: %f [pm^2]\n', area(shp));

subplot(3,2,4);
s = scatter(r12, theta, 'o');
s.MarkerEdgeColor = [0.000, 0.447, 0.741];
s.MarkerEdgeAlpha = 0.6;
xlabel('r_{CO} (pm)');
ylabel('\theta (deg)');
axis tight;
box on;

hold on;
[k,v] = convhull(r12, theta);
p = plot(r12(k), theta(k), 'r:', 'LineWidth', 1.5);
p.Color = [0.635, 0.078, 0.184, 0.6];
fprintf('r12*theta convex hull area: %f [pm*deg]\n', v);

hold on;
shp = alphaShape(r12, theta, 22);
bf = boundaryFacets(shp);
p = plot(shp.Points(bf,1), shp.Points(bf,2), '-', 'LineWidth', 1.5);
p.Color = [0.494, 0.184, 0.556, 0.6];
fprintf('r12*theta alpha shape area: %f [pm*deg]\n', area(shp));

subplot(3,2,6);
s = scatter(r23, theta, 'o');
s.MarkerEdgeColor = [0.000, 0.447, 0.741];
s.MarkerEdgeAlpha = 0.6;
xlabel('r_{CS} (pm)');
ylabel('\theta (deg)');
axis tight;
box on;

hold on;
[k,v] = convhull(r23, theta);
p = plot(r23(k), theta(k), 'r:', 'LineWidth', 1.5);
p.Color = [0.635, 0.078, 0.184, 0.6];
fprintf('r23*theta convex hull area: %f [pm*deg]\n', v);

hold on;
shp = alphaShape(r23, theta, 30);
bf = boundaryFacets(shp);
p = plot(shp.Points(bf,1), shp.Points(bf,2), '-', 'LineWidth', 1.5);
p.Color = [0.494, 0.184, 0.556, 0.6];
fprintf('r23*theta alpha shape area: %f [pm*deg]\n', area(shp));

[k,v] = convhull(r12, r23, theta);
shp = alphaShape(r12, r23, theta);
fprintf('r12*r23*theta convex hull volume: %f [pm^2*deg]\n', v);
fprintf('r12*r23*theta alpha shape volume: %f [pm^2*deg]\n', volume(shp));

end
