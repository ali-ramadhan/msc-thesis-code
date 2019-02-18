function exploration1(error)

% Corresponds to a geometry of
% r_CO=131.677 pm, r_CS=187.161 pm, theta=169.42 deg 
% (bond lengths stretched by ~0.3 A, slightly bent. This is roughly in the
% center of physically reasonable geometries in phase space.)
p = [5.4289900e-22 -3.5614700e-22 3.7813800e-23 1.8521600e-22 ...
     1.7229300e-22 8.4222500e-23 -7.0915100e-22 1.6204400e-22 ...
     -1.1635100e-22];

pAll = generateAllPlusMinusMomentumVectors(p, error);
 
masses = [16 12.011 32.06];
charges = [2 2 2];
filenamePrefix = strcat('EXP1_', num2str(error));

multiStartTriatomic(pAll, masses, charges, filenamePrefix, 1, 1);
end