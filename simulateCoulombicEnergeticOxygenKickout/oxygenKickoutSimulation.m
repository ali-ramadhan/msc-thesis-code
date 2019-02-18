function out = oxygenKickoutSimulation(E_mu, E_sigma, n)

eV = 1.602e-19;
E = eV*normrnd(E_mu, E_sigma, [n 1]);
E = E - 2*min(0, E);

r_12  = normrnd(120, 5, [n 1]);
r_23  = normrnd(160, 7.55, [n 1]);
theta = normrnd(170, 5, [n 1]);

r_12 = r_12 - 2*min(0,r_12-110);
r_23 = r_23 - 2*min(0,r_23-150);
theta = theta - 2*max(0, theta-180);

r_12 = 1e-12*r_12;
r_23 = 1e-12*r_23;

out = zeros(n,9);

for i=1:1:n
    if rem(i, 1000) == 0
        fprintf('%d/%d\n', i, n);
    end
    
    out(i,:) = simulateCoulombicEnergeticOxygenKickout(r_12(i), r_23(i), theta(i), E(i));
end

end