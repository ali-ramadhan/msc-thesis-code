function out = generateAllPlusMinusMomentumVectors(p, error)
  lowerBound = 1-error;
  upperBound = 1+error;
  pm = [lowerBound; upperBound];
  n = size(pm,1);

  [i1,i2,i3,i4,i5,i6,i7,i8,i9] = ndgrid(1:n, 1:n, 1:n, 1:n, 1:n, 1:n, 1:n, 1:n, 1:n);
  prod = [pm(i1,:), pm(i2,:), pm(i3,:), pm(i4,:), pm(i5,:), pm(i6,:), pm(i7,:), pm(i8,:), pm(i9,:)];

  nOut = size(prod,1);

  out = prod .* (ones(nOut, 1) * p);
end
