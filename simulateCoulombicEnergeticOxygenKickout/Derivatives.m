function out = Derivatives (Time, Par)
% Derivatives.m (c) J.P. Brichta 2007 - 2009
% Inputs:
% --- Time : a 1x2 row vector containing [InitialTime, FinalTime]
% --- Par : a 1x18 row vector containing position and momemtum parameters of the three particles.
%           The particles' co-ordinate components are given as [x1, y1, z1, x2, y2, z2, x3, y3, z3, px1, ..., pz3].
% Output:
% --- out : a nx19 array where each row contains [Time, Position[1x9], Momentum[1x9]] in the same format as Par. In
%           practice, only the final row is utilized to evaluate the final conditions of the system.
% Notes: All units are SI and are stated for clarity at the end of variable assignment statements.

% Constants.
amu = 1.66053886e-27; % [kg], 1 atomic mass unit
e   = 1.60217646e-19; % [C], 1 elementary charge
k   = 8.987551e9;     % [N m^2 C^-2], Electrostatic constant

% System specific values
% Only change the "natural" number. For instance, if you want to change the charge on particle 3 to be 2+, the
% assignment should be q3 = 2*e. In principle, negative charge states can be declared, but this may lead to unusual
% results in the case of mixed (postive and negative) ions.

m1 = 15.9994*amu; % [kg], Mass of Oxygen
m2 = 12*amu;  % [kg], Mass of Carbon
m3 = 32.065*amu;  % [kg], Mass of Sulphur

q1 = 2*e; % [C], Charge on Oxygen
q2 = 1*e; % [C], Charge on Carbon
q3 = 1*e; % [C], Charge on Sulphur

% Calculate the distance between ions. Note that this quantity does not preserve vector direction.
r12 = ((Par(1)-Par(4))^2 + (Par(2)-Par(5))^2 + (Par(3)-Par(6))^2)^0.5; % [m]
r13 = ((Par(1)-Par(7))^2 + (Par(2)-Par(8))^2 + (Par(3)-Par(9))^2)^0.5; % [m]
r23 = ((Par(4)-Par(7))^2 + (Par(5)-Par(8))^2 + (Par(6)-Par(9))^2)^0.5; % [m]

% ParDot is a column vector with components [vx1; vy1; vz1; vx2; ... vz3; p'x1; ... p'z3]. These quantities are produced
% by taking the first derivative of the Hamiltonian with respect to the appropriate variable. For a complete derivation,
% refer to Brichta et al, Computer Physics Communications, vol. 180 (2009) 197-200, equations 1-7.

ParDot = [Par(10)./m1; Par(11)./m1; Par(12)./m1; ...
          Par(13)./m2; Par(14)./m2; Par(15)./m2; ...
          Par(16)./m3; Par(17)./m3; Par(18)./m3; ...
      
          k*q1*q2*(Par(1)-Par(4))/r12^3 + k*q1*q3*(Par(1)-Par(7))/r13^3; ...
          k*q1*q2*(Par(2)-Par(5))/r12^3 + k*q1*q3*(Par(2)-Par(8))/r13^3; ...
		  k*q1*q2*(Par(3)-Par(6))/r12^3 + k*q1*q3*(Par(3)-Par(9))/r13^3; ...
      
          k*q2*q1*(Par(4)-Par(1))/r12^3 + k*q2*q3*(Par(4)-Par(7))/r23^3; ...
		  k*q2*q1*(Par(5)-Par(2))/r12^3 + k*q2*q3*(Par(5)-Par(8))/r23^3; ...
		  k*q2*q1*(Par(6)-Par(3))/r12^3 + k*q2*q3*(Par(6)-Par(9))/r23^3; ...
      
          k*q3*q1*(Par(7)-Par(1))/r13^3 + k*q3*q2*(Par(7)-Par(4))/r23^3; ...
		  k*q3*q1*(Par(8)-Par(2))/r13^3 + k*q3*q2*(Par(8)-Par(5))/r23^3; ...
          k*q3*q1*(Par(9)-Par(3))/r13^3 + k*q3*q2*(Par(9)-Par(6))/r23^3];

out = ParDot;