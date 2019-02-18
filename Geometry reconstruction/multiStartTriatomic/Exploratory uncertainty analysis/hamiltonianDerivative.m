function out = hamiltonianDerivative(time, par)
% Inputs:
% * time: a 1x2 row vector containing [InitialTime, FinalTime]
% * par : a 1x18 row vector containing position and momemtum parameters of the
%         three particles. The particles' co-ordinate components are given as
%         [x1, y1, z1, x2, y2, z2, x3, y3, z3, px1, ..., pz3].
% Output:
% * out : an nx19 array where each row contains [Time, Position[1x9], Momentum[1x9]]
%         in the same format as par. In practice, only the final row is utilized
%         to evaluate the final conditions of the system.
% Notes: All units are SI.

% Constants.
amu = 1.66053886e-27; % [kg], 1 atomic mass unit
e   = 1.60217646e-19; % [C], 1 elementary charge
k   = 8.987551e9;     % [N m^2 C^-2], electrostatic constant

% Masses and charges
m1 = amu*par(19); m2 = amu*par(20); m3 = amu*par(21);
q1 = e*par(22);   q2 = e*par(23);   q3 = e*par(24);

% Calculate the distance between ions. Note that this quantity does not preserve
% vector direction.
r12 = ((par(1)-par(4))^2 + (par(2)-par(5))^2 + (par(3)-par(6))^2)^0.5; % [m]
r13 = ((par(1)-par(7))^2 + (par(2)-par(8))^2 + (par(3)-par(9))^2)^0.5; % [m]
r23 = ((par(4)-par(7))^2 + (par(5)-par(8))^2 + (par(6)-par(9))^2)^0.5; % [m]

% parDot is a column vector with components [vx1; vy1; vz1; vx2; ... vz3;
% p'x1; ... p'z3]. These quantities are produced by taking the first derivative
% of the Hamiltonian with respect to the appropriate variable. For a complete
% derivation, refer to Brichta et al, Computer Physics Communications, vol. 180
% (2009) 197-200, equations 1-7.

parDot = [par(10)./m1; par(11)./m1; par(12)./m1; ...
          par(13)./m2; par(14)./m2; par(15)./m2; ...
          par(16)./m3; par(17)./m3; par(18)./m3; ...

          k*q1*q2*(par(1)-par(4))/r12^3 + k*q1*q3*(par(1)-par(7))/r13^3; ...
          k*q1*q2*(par(2)-par(5))/r12^3 + k*q1*q3*(par(2)-par(8))/r13^3; ...
		  k*q1*q2*(par(3)-par(6))/r12^3 + k*q1*q3*(par(3)-par(9))/r13^3; ...

          k*q2*q1*(par(4)-par(1))/r12^3 + k*q2*q3*(par(4)-par(7))/r23^3; ...
		  k*q2*q1*(par(5)-par(2))/r12^3 + k*q2*q3*(par(5)-par(8))/r23^3; ...
          k*q2*q1*(par(6)-par(3))/r12^3 + k*q2*q3*(par(6)-par(9))/r23^3; ...

          k*q3*q1*(par(7)-par(1))/r13^3 + k*q3*q2*(par(7)-par(4))/r23^3; ...
          k*q3*q1*(par(8)-par(2))/r13^3 + k*q3*q2*(par(8)-par(5))/r23^3; ...
          k*q3*q1*(par(9)-par(3))/r13^3 + k*q3*q2*(par(9)-par(6))/r23^3];

out = [parDot; par(19:24)];
end
