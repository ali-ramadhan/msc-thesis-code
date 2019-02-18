function out = stepwiseSimulation(E1norm, rotationPeriod, t2, E2norm, tEnd, delta_t, n)

eV = 1.602e-19;

CSBondEnergy = E1norm(1);
DeltaCSBondEnergy = E1norm(2);
COBondEnergy = E2norm(1);
DeltaCOBondEnergy = E2norm(2);

E1 = eV*normrnd(CSBondEnergy, DeltaCSBondEnergy, [n 1]);
E2 = eV*normrnd(COBondEnergy, DeltaCOBondEnergy, [n 1]);

E1 = E1 - 2*min(0, E1);
E2 = E2 - 2*min(0, E2);

out = zeros(n,9);

for i=1 : 1 : n
    if rem(i, 100) == 0
        fprintf('%d/%d\n', i, n);
    end
    
    out(i,:) = simulateSequentialBondBreaking(E1(i), rotationPeriod, t2, E2(i), tEnd, delta_t);
end

end