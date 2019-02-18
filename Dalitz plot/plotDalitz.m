function [X,Y] = plotDalitz(momentumFileName)
momenta = load(momentumFileName);
momenta = momenta*1e-22;

for i = 1:size(momenta,1)
    momenta(i,:) = removeCOMMotion(momenta(i,:));
    momenta(i,:) = rotateMomentum2(momenta(i,:));
end

p_O = momenta(:,1:3);
p_C = momenta(:,4:6);
p_S = momenta(:,7:9);

amu = 1.66053886e-27;
e   = 1.602176565e-19;

m_O = 15.999*amu;
m_C = 12.011*amu;
m_S = 32.065*amu;

E_O = sum(p_O.^2, 2) / (2*m_O);
E_C = sum(p_C.^2, 2) / (2*m_C);
E_S = sum(p_S.^2, 2) / (2*m_S);
E_total = E_O + E_C + E_S;

hist(E_O,50)
figure;
hist(E_C,50)
figure;
hist(E_S,50)

E = [E_O E_C E_S E_total];
E = (1/e)*E;

% size(E)
% mean(E)
% 
% for r = 1:size(E,1)
%     if E(r,4) < 30 || E(r,4) > 200
%         E(r,:) = [0 0 0 0];
%     end
% end
% E(~any(E,2), :) = [];
% 
% size(E)
% mean(E)
% min(E)
% max(E)

E_O = E(:,1);
E_C = E(:,2);
E_S = E(:,3);
E_total = E(:,4);

epsilon_O = E_O ./ E_total;
epsilon_C = E_C ./ E_total;
epsilon_S = E_S ./ E_total;

X = (epsilon_O - epsilon_S) / sqrt(3);
Y = epsilon_C - (1/3);

save 'momentaAr6+X_rotated.txt' X -ascii -double
save 'momentaAr6+Y_rotated.txt' Y -ascii -double
end

%-----------------------------------------------------------------------------------------------------------------------
% removeCOMMOtion takes a momentum triple, [p_1 p_2 p_3] for example, and returns the triple in the same order except
% with the center of mass motion removed.
function out = removeCOMMotion(momentum)
% We split our momentum triple into X,Y,Z components.
p_X = momentum(1:3:7);
p_Y = momentum(2:3:8);
p_Z = momentum(3:3:9);

mass = [12.0107 15.9994 32.065]; % OCS
massSum = sum(mass);

% We eliminate COM motion to put the particles in the COM frame of motion.
p_X = p_X - sum(p_X) .* (mass)/massSum;
p_Y = p_Y - sum(p_Y) .* (mass)/massSum;
p_Z = p_Z - sum(p_Z) .* (mass)/massSum;

% Putting the vectors back into the original format.
momentum(1:9) = reshape([p_X; p_Y; p_Z], 1, 9);
out = momentum;
end
%-----------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------
% rotateMomentum2 rotates the momentum to eliminate the z-component. We change format 123->213 (e.g. OCS->COS) because
% rotateMomentum expects 213/COS.
function out = rotateMomentum2(momentum)
p_1 = momentum(1:3);
p_2 = momentum(4:6);
p_3 = momentum(7:9);

momentum = [p_1 p_2 p_3];
momentum = rotateMomentum(momentum);
momentum = [momentum(4:6) momentum(1:3) momentum(7:9)];

out = momentum(1:9);
end
%-----------------------------------------------------------------------------------------------------------------------