h = 0.23;
E = 1900;

amu = 1.66053886e-27;
e = 1.60217646e-19;

m = linspace(1,100);
t1 = 1e9 * sqrt( (2*h) .* (amu*m) / (e*E));
t2 = 1e9 * sqrt( (2*h) .* (amu*m) / (2*e*E));

plot(m, t1, '-r', m, t2, '-g', 'LineSmoothing', 'on');
legend('q = 1', 'q = 2', 'Location', 'SouthEast');
title('Time of Flight');
xlabel('Mass (amu)');
ylabel('Time (ns)');
set(gca,'Color',[0.39 0.47 0.64]);
grid on;