function out = simulateCoulombicEnergeticOxygenKickout(r_12, r_23, theta, E_O_kickout)
amu  = 1.66053886e-27;
e    = 1.60217646e-19;
k_e  = 8.9875517873681764e9;

m_O = 16*amu;
m_C = 12*amu;
m_S = 32*amu;

r_O = [-r_12, 0];
r_C = [0, 0];
r_S = [r_23*cosd(180-theta), r_23*sind(180-theta)];
r_CS = (m_S*r_S + m_C*r_C) / (m_S + m_C); % Using COM.

O_speed = sqrt(2*E_O_kickout / m_O);
v_O  = O_speed * ((r_O - r_CS) / norm(r_O - r_CS)); % Kick O out in direction away from CS COM.

OCS   = [r_O 0 r_C 0 r_S 0];
p_OCS = [m_O*v_O 0, 0 0 0, 0 0 0];
[t,y] = IonVelocities([OCS p_OCS]);

p_O = y(end,10:12);
p_C = y(end,13:15);
p_S = y(end,16:18);

% R_O = y(:,1:2);
% R_C = y(:,4:5);
% R_S = y(:,7:8);
% 
% plot(R_O(:,1), R_O(:,2), '-ro', R_C(:,1), R_C(:,2), '-ko', R_S(:,1), R_S(:,2), '-yo', 'LineSmoothing', 'on');
% legend('O', 'C', 'S', 'Location', 'SouthEast');
% axis([-10e-10 10e-10 -10e-10 10e-10]);
% title('Ion position (short timescale)');
% xlabel('X (m)');
% ylabel('Y (m)');
% set(gca,'Color',[0.39 0.47 0.64]);
% grid on;

out = [p_O p_C p_S];
end

% IonVelocities: Calculate velocities of ion given initial parameters
% Usage: IonVelocities(InitialConditions) where InitialConditions is a vector 
% of the form [x1, y1, z1, x2, ..., z3].
function [t,y] = IonVelocities (InitialConditions)
    options = odeset('AbsTol', 1e-27, 'RelTol', 1e-6, 'InitialStep', 1e-18);
    [t,y] = ode45('Derivatives', [0 1e-11], InitialConditions, options);
end