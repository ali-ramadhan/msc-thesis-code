function out = rotateMomentum (momenta)
% Input:
% * momenta: a 1x9 array of p_x, p_y, p_z triples.
%   By convention, the first triple set is the center atom's momentum vectors.
% Output:
% * out: a 1x11 arrary of p_x, p_y, 0 triples, theta_v, chi.
%   By convention, the first triple set is the oxygen 1 momentum vectors.
%
% Rotate asymptotic momentum vectors of Coulomb explosion such that output
% vectors lie in the x-y plane.  Of course, this transformation preserves
% geometries.
%
% The second part of the transformation is more for neatness.  With this
% transformation we rotate everything in the x-y plane so that the x-axis
% bisects the angle between the two oxygen ions. The carbon just comes along for
% the ride.

% We want to put the momentum in the form [p_x1, p_y1, p_z1; p_x2, p_y2, p_z2; p_x3, p_y3, p_z3].
momenta = reshape(momenta, 3, 3)';

% If all the vectors are already in the xy-plane, just return the vectors.
z_components = momenta(:,3);
if z_components == [0; 0; 0]
  theta_v = acos(dot(momenta(2,:), momenta(3,:)) / norm(momenta(2,:)) / norm(momenta(3,:))); % WTF: two divisions?
  chi = acos(dot(momenta(1,:), momenta(2,:)-momenta(3,:))/norm(momenta(1,:))/norm(momenta(2,:)-momenta(3,:)));
  out =  [momenta(2,:), momenta(1,:), momenta(3,:), theta_v, chi];
  return
end

% First, check that all three momentum vectors form a plane.
% If this is the case, the determinant of the momentum vectors will be zero (or
% a very small number).

if det(momenta) < 1e-50
  % The normal vector of the plane is the cross product of two of the vectors.
  normal = cross(momenta(1,:), momenta(2,:));

  % Normalise the normal vector.
  normal = normal / norm(normal);

  % The components of the normal vector (A,B,C) define the plane.  By
  % construction, the plane goes through the origin.  We must now find the angle
  % this plane makes with the x-y plane, which is defined by the normal vector
  % normal_xy = [0,0,1].  This angle is the dihedral angle.
  normal_xy = [0 0 1];
  dihedral = acos(dot(normal, normal_xy)); % [rad]

  % The two planes intersect in a line defined by crossing the two normal
  % vectors.
  intersection = cross(normal, normal_xy);

  %if intersection - zeros(1,3) ~= 0 % this is only the case if
  intersection = intersection / norm(intersection);

  % Insert all of the relevant information into the matrix which performs the
  % plane rotation.
  momenta = rotatePlane(momenta, dihedral, intersection);

  % Determine theta_v, the angle between the sulpher(??) fragments
  theta_v = acos(dot(momenta(2,:), momenta(3,:))/norm(momenta(2,:))/norm(momenta(3,:)));

  % Determine phi, the minimum angle between one of the oxygens and the x-axis.
  % phi = abs(atan2(momenta(2,2), momenta(2,1)));

  phi = abs(atan(momenta(2,2)/momenta(2,1)));

  if momenta(2,1) >= 0 && momenta(2,2) >= 0 % first quadrant
      phi = phi;
  elseif momenta(2,1) < 0 && momenta(2,2) >= 0 % second quadrant
      phi = pi - phi;
  elseif momenta(2,1) < 0 && momenta(2,2) < 0 % third quadrant
      phi = pi + phi;
  else % fourth quadrant
      phi = 2*pi - phi;
  end

  % Rotate everything (clockwise) through the angle phi
  M = [cos(-phi), -sin(-phi), 0; sin(-phi), cos(-phi), 0; 0, 0, 1];
  for i = 1:3
    momenta(i,:) = (M*(momenta(i,:)'))';
  end

  % Flip in the y-axis (if necessary) such that the second end atom sits in the
  % +x half plane.
  if momenta(3,2) < 0
    momenta(1:2:3,2) = -momenta(1:2:3,2);
  end

  if rand(1) > 0.5
      chi = acos(dot(momenta(1,:), momenta(2,:)-momenta(3,:))/norm(momenta(1,:))/norm(momenta(2,:)-momenta(3,:)));
  else
      chi = acos(dot(momenta(1,:), momenta(3,:)-momenta(2,:))/norm(momenta(1,:))/norm(momenta(3,:)-momenta(2,:)));
  end

  out = [momenta(2,:), momenta(1,:), momenta(3,:), theta_v, chi];

  else
    out = zeros(1,11);
  end
end

% Rotate the momentum vectors which define the plane Ax + By + Cz = 0 by the
% angle dihedral in the line of intersection defined by the direction cosines a,
% b, and c made with the plane z=0. This makes the momentum vectors exist only
% in the x-y plane while retaining their configuration and magnitude. The
% dihedral angle is in radians.
function out = rotatePlane (momenta, dihedral, intersection)
  lineMagnitude = norm(intersection);
  a = intersection(1) / lineMagnitude;
  b = intersection(2) / lineMagnitude;
  c = intersection(3) / lineMagnitude;

  M = [a^2*(1-cos(dihedral))+cos(dihedral),   a*b*(1-cos(dihedral))-c*sin(dihedral), a*c*(1-cos(dihedral))+b*sin(dihedral); ...
       a*b*(1-cos(dihedral))+c*sin(dihedral), b^2*(1-cos(dihedral))+cos(dihedral),   b*c*(1-cos(dihedral))-a*sin(dihedral); ...
       a*c*(1-cos(dihedral))-b*sin(dihedral), b*c*(1-cos(dihedral))+a*sin(dihedral), c^2*(1-cos(dihedral))+cos(dihedral)];

  for i = 1:3
    newVector(i,:) = (M*(momenta(i,:)'))';
  end

  out = newVector;
end
