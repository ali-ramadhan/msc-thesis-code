% simulateMomenta takes in
% * geometries: a matrix where each row is of the form [r_12 r_23 theta].
%   r_12 and r_23 should be given in SI units [m] and theta in [deg].
% * masses:     a vector [m1 m2 m3] with the atomic masses in amu.
% * charges:    a vector [q1 q2 q3] with the atomic charges in units of e.
function out = simulateMomentum(geometry, masses, charges)
  out = zeros(1,12);

  r_12  = geometry(1);
  r_23  = geometry(2);
  theta = geometry(3);

  % Place the first atom to the left of central atom.
  x_1 = -r_12;
  y_1 = 0;
  z_1 = 0;

  % Place the central atom at the origin.
  x_2 = 0;
  y_2 = 0;
  z_2 = 0;

  % Place the third atom to the right of the central taking into the account
  % the angle between the two bond lengths.
  x_3 = r_23 * cosd(180 - theta);
  y_3 = r_23 * sind(180 - theta);
  z_3 = 0;

  % Calculate the momentum that such a geometry would produce.
  g = [x_1 y_1 z_1 x_2 y_2 z_2 x_3 y_3 z_3];
  p_0 = zeros(1,9);
  p = ionVelocities([g p_0], masses, charges);
  out = [r_12 r_23 theta p(1:3) p(4:6) p(7:9)];

  % Extract each atom's momentum into 2D vectors in preparation to
  % rotate.
  p_1 = out(4:5);
  p_2 = out(7:8);
  p_3 = out(10:11);

  % Put each momentum into column vector form so we can use matrix
  % multiplication.
  p_1 = p_1'; p_2 = p_2'; p_3 = p_3';

  % Calculate the angle between the central atom and the +x-axis then
  % rotate the three momentum vectors back towards the origin by that
  % much so that the central's momentum is always along the +x-axis.
  theta_2x = atan2(p_2(2), p_2(1));
  R = [cos(-theta_2x) -sin(-theta_2x); sin(-theta_2x) cos(-theta_2x);];
  p_1 = R*p_1;
  p_2 = R*p_2;
  p_3 = R*p_3;

  % Put everything back into a row vector.
  p_1 = p_1'; p_2 = p_2'; p_3 = p_3';

  % Set the z components (and also y in case of carbon) to zero so
  % they all have exactly the same value rather than 0.0000 and
  % -0.0000, etc.
  p_1 = [p_1 0];
  p_2 = [p_2(1) 0 0];
  p_3 = [p_3 0];

  out = [r_12 r_23 theta p_1 p_2 p_3];
end

% IonVelocities: Calculate velocities of ion given initial parameters
% Usage: IonVelocities(InitialConditions) where InitialConditions is a vector
% of the form [x1, y1, z1, x2, ..., z3].
function out = ionVelocities (initialConditions, masses, charges)
  options = odeset('AbsTol', 1e-27, 'RelTol', 1e-6, 'InitialStep', 1e-18);
  [t,y] = ode45('hamiltonianDerivative', [0 1e-11], [initialConditions masses charges], options);
  dlmwrite('myFile.txt',t)
  out = y(size(t,1), 10:18);
end
