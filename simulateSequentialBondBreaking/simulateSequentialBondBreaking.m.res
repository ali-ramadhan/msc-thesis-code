function out = simulateSequentialBondBreaking(amu_1, amu_2, amu_3, q_1, q_2, q_3, r_12, r_23, theta, t_secondBondBreaks, delta_t)

amu = 1.66053886e-27;
e = 1.60217646e-19;
k_e  = 8.9875517873681764e9;

m_1 = amu_1 * amu;
m_2 = amu_2 * amu;
m_3 = amu_3 * amu;

q_1 = q_1 * e;
q_2 = q_2 * e;
q_3 = q_3 * e;

r_1 = [-r_12, 0, 0];
r_2 = [0, 0, 0];
r_3 = [r_23*cosd(180-theta), r_23*sind(180-theta), 0];

R_1 = [r_1];
R_2 = [r_2];
R_3 = [r_3];

v_1 = 0;
v_2 = 0;
v_3 = 0;

V_1 = [v_1];
V_2 = [v_2];
V_3 = [v_3];

t = 0;

% Time evolve the system under a Coulombic force with ions 2 and 3 bonded (constrained to have their positions
% constant relative to each other) but separated from ion 1. We account for torque and force.
while t < t_secondBondBreaks
    a_1 = ((k_e*q_1) / m_1) * ( (q_2*(r_2-r_1))/(norm(r_2-r_1)^3) + (q_3*(r_3-r_1))/(norm(r_3-r_1)^3));
    a_2 = ((k_e*q_2) / m_2) * ( (q_1*(r_1-r_2))/(norm(r_1-r_2)^3) + (q_3*(r_3-r_2))/(norm(r_3-r_2)^3));
    a_3 = ((k_e*q_3) / m_3) * ( (q_1*(r_1-r_3))/(norm(r_1-r_3)^3) + (q_2*(r_2-r_3))/(norm(r_2-r_3)^3));
    
    v_1 = v_1 + a_1*delta_t;
    v_2 = v_2 + a_2*delta_t;
    v_3 = v_3 + a_3*delta_t;
    
    V_1(end+1) = v_1;
    V_2(end+1) = v_2;
    V_3(end+1) = v_3;
    
    r_1 = r_1 + v_1*delta_t;
    r_2 = r_2 + (v_2 + v_3)*delta_t;
    r_3 = r_3 + (v_2 + v_3)*delta_t;
    
    R_1(end+1) = r_1;
    R_2(end+1) = r_2;
    R_3(end+1) = r_3;  
    
    t = t + delta_t;
end

end 