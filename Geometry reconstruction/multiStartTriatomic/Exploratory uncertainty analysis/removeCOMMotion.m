% removeCOMMOtion takes a momentum triple, [p_1 p_2 p_3] and returns the triple
% in the same order with the center of mass motion removed.
% masses is the mass of each of the atoms (O, then C, then S) in amu.
function out = removeCOMMotion(momentum, masses)
  % We split our momentum triple into X,Y,Z components.
  p_X = momentum(1:3:7);
  p_Y = momentum(2:3:8);
  p_Z = momentum(3:3:9);

  masses = [15.9994 12.0107 32.065];  % OCS
  massSum = sum(masses);

  % We eliminate COM motion to put the particles in the COM frame of motion.
  p_X = p_X - sum(p_X) .* (masses)/massSum;
  p_Y = p_Y - sum(p_Y) .* (masses)/massSum;
  p_Z = p_Z - sum(p_Z) .* (masses)/massSum;

  % Putting the vectors back into the original format.
  momentum(1:9) = reshape([p_X; p_Y; p_Z], 1, 9);
  out = momentum;
end
