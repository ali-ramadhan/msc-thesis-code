function out = simulateSequentialBondBreaking(E_firstBondBreaks, rotationPeriod, t_secondBondBreaks, E_secondBondBreaks, t_endNG, delta_t)

global t_end
t_end = t_endNG;

amu = 1.66053886e-27;
e = 1.60217646e-19;
k_e  = 8.9875517873681764e9;

m_O = 16*amu;
m_C = 200*amu;
m_S = 32*amu;
m_CO = m_C + m_O;

% Set up the molecule in the initial ground state geometry.
r_12  = 115.78e-12;
r_23  = 156.01e-12;
theta = 170;

r_O = [-r_12, 0];
r_C = [0, 0];
r_S = [r_23*cosd(180-theta), r_23*sind(180-theta)];

r_CO = (m_O*r_O + m_C*r_C) / (m_O + m_C); % Using COM.
r_CO2C = r_C - r_CO;
r_CO2O = r_O - r_CO;

v_CO = [0,0];
v_S  = [0,0];

% Break the first bond (the CS bond) and deposit the bond dissociation energy into the system as kinetic energy,
% distributing it such that momentum is conserved.
CO_speed = sqrt( 2 * ( m_S / (m_CO*m_S + m_CO^2) ) * E_firstBondBreaks );
S_speed  = sqrt( 2 * ( m_CO / (m_CO*m_S + m_S^2) ) * E_firstBondBreaks );

v_CO = v_CO + CO_speed * ((r_CO - r_S) / norm(r_CO - r_S));
v_S  = v_S + S_speed * ((r_S - r_CO) / norm(r_S - r_CO));

% Start at zero time.
t = 0;
i = 1;

% Ionize each fragment once to start this whole thing off.
q_CO = 2*e;
q_S  = e;

% Record the positions and velocities of each fragment for all time. These storage vectors use capital letters.
n = ceil(t_secondBondBreaks/delta_t);
R_C  = zeros(n,2); R_C(i,:)  = r_C;
R_O  = zeros(n,2); R_O(i,:)  = r_O;
R_CO = zeros(n,2); R_CO(i,:) = r_CO;
R_S  = zeros(n,2); R_S(i,:)  = r_S;

V_CO = zeros(n,2); V_CO(i,:) = v_CO;
V_S  = zeros(n,2); V_S(i,:)  = v_S;

T = [t];

omega = 2*pi / rotationPeriod; % [rad]

% Time evolve the system under a Coulombic force with the CO fragment still bonded. The fragment rotates at some contant
% angular velocity about its COM.
while t < t_secondBondBreaks
    % Calculate the rotational motion of the CO fragment and rotate it a little bit. This doesn't change the position
    % of the CO fragment.
%     r_C_COM = r_C - r_CO;
%     r_O_COM = r_O - r_CO;
%     
%     polar_r_C     = norm(r_C_COM);
%     polar_theta_C = atan2(r_C_COM(2), r_C_COM(1));
%     polar_r_O     = norm(r_O_COM);
%     polar_theta_O = atan2(r_O_COM(2), r_O_COM(1));
%     
%     polar_theta_C = polar_theta_C + omega*delta_t;
%     polar_theta_O = polar_theta_O + omega*delta_t;
%     
%     r_C_COM = polar_r_C * [cos(polar_theta_C), sin(polar_theta_C)];
%     r_O_COM = polar_r_O * [cos(polar_theta_O), sin(polar_theta_O)];
%     
%     r_C = r_CO + r_C_COM;
%     r_O = r_CO + r_O_COM;
    
    % Calculate fragment acceleration due to the Coulomb repulsion and move the fragments a little bit.
    a_CO = -((k_e*q_CO) / m_CO) * ( (q_S*(r_S-r_CO)) / (norm(r_S-r_CO)^3) );
    a_S  = -((k_e*q_S) / m_S) * ( (q_CO*(r_CO-r_S)) / (norm(r_CO-r_S)^3) );
    
    v_CO = v_CO + a_CO*delta_t;
    v_S  = v_S  + a_S*delta_t;
    
    r_CO = r_CO + v_CO*delta_t;
    r_O  = r_O + v_CO*delta_t;
    r_C  = r_C + v_CO*delta_t;
    r_S  = r_S  + v_S*delta_t;
    
    % Record their new positions.
	i = i+1;
    R_C(i,:)  = r_CO + r_CO2C;
	R_O(i,:)  = r_CO + r_CO2O;
	R_CO(i,:) = r_CO;
	R_S(i,:)  = r_S;

	V_CO(i,:) = v_CO;
	V_S(i,:)  = v_S;
    
    t = t + delta_t;
    T(end+1) = t;
    
%     if rem(i, 1000) == 0
%         fprintf('%d/%d\n', t, t_secondBondBreaks);
%     end
end

% Break the second bond (the CO bond) and deposit the bond dissociation energy into the now separated rigid rotor as
% kinetic energy directed away from the CO COM position such that the fragments explode away from each other. It's
% distributed such that momentum is conserved. This kinectic energy is added to the kinetic energy due to COM
% translational motion and the rotation of each atom about the COM.
C_speedFromBondBreaking = sqrt( 2 * ( m_O / (m_O*m_C + m_C^2) ) * E_secondBondBreaks );
O_speedFromBondBreaking = sqrt( 2 * ( m_C / (m_O*m_C + m_O^2) ) * E_secondBondBreaks );

% Rotate the CO rotor by a random angle [0,2*pi] to simulate random explosion orientation.
theta_rand = 2*pi*rand(1);
R_theta = [cos(theta_rand) -sin(theta_rand); sin(theta_rand) cos(theta_rand);]; % 2D xy rotation matrix.

r_C_COM = (r_C - r_CO)*R_theta;
r_O_COM = (r_O - r_CO)*R_theta;

r_C = r_CO + r_C_COM;
r_O = r_CO + r_O_COM;

v_C = v_CO + C_speedFromBondBreaking * (r_C_COM / norm(r_C_COM));
v_O = v_CO + O_speedFromBondBreaking * (r_O_COM / norm(r_O_COM));

omegaVector = [0, 0, omega]; % omega vector with right hand rule convention
v_C_fromRotation = cross(omegaVector, [r_C_COM 0]);
v_O_fromRotation = cross(omegaVector, [r_O_COM 0]);

v_C = v_C + v_C_fromRotation(1:2);
v_O = v_O + v_O_fromRotation(1:2);

R_O = [R_O; r_O];
R_C = [R_C; r_C];
R_S = [R_S; r_S];
    
V_CO = [V_CO; v_CO];
V_S  = [V_S; v_S];

V_O = [v_O];
V_C = [v_C];

OCS   = [r_O 0 r_C 0 r_S 0];
p_OCS = [m_O*v_O 0, m_C*v_C 0, m_S*v_S 0];
[t,y] = IonVelocities([OCS p_OCS]);

R_O = [R_O; y(:,1:2)];
R_C = [R_C; y(:,4:5)];
R_S = [R_S; y(:,7:8)];

p_O = y(end,10:12);
p_C = y(end,13:15);
p_S = y(end,16:18);

out = [p_O p_C p_S];

% subplot(1,2,1);
% plot(R_O(:,1), R_O(:,2), '-ro', R_C(:,1), R_C(:,2), '-ko', R_CO(:,1), R_CO(:,2), '-bo', R_S(:,1), R_S(:,2), '-yo', 'LineSmoothing', 'on');
% legend('O', 'C', 'CO', 'S', 'Location', 'SouthEast');
% axis([-5e-10 5e-10 -5e-10 5e-10]);
% title('Ion position (short timescale)');
% xlabel('X (m)');
% ylabel('Y (m)');
% set(gca,'Color',[0.39 0.47 0.64]);
% grid on;
% 
% subplot(1,2,2);
% plot(R_O(:,1), R_O(:,2), '-ro', R_C(:,1), R_C(:,2), '-ko', R_CO(:,1), R_CO(:,2), '-bo', R_S(:,1), R_S(:,2), '-yo', 'LineSmoothing', 'on');
% legend('O', 'C', 'CO', 'S', 'Location', 'SouthEast');
% title('Ion position (long timescale)');
% xlabel('X (m)');
% ylabel('Y (m)');
% set(gca,'Color',[0.39 0.47 0.64]);
% grid on;
end 

% IonVelocities: Calculate velocities of ion given initial parameters
% Usage: IonVelocities(InitialConditions) where InitialConditions is a vector 
% of the form [x1, y1, z1, x2, ..., z3].
function [t,y] = IonVelocities (InitialConditions)
    global t_end
    options = odeset('AbsTol', 1e-27, 'RelTol', 1e-6, 'InitialStep', 1e-18);
    [t,y] = ode45('Derivatives', [0 t_end], InitialConditions, options);
end