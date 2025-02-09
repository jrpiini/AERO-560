%% AERO 560 HW 3 - Joseph Piini

clear 
close all
clc

% addpath('C:\Users\josep\OneDrive - Cal Poly\Desktop\AERO 560\Functions')

%% Problem 3

% Orbit Parameters

mu = 398600.4418;

h = 55759; % km^2/s
ecc = 0.001;
RAAN = 10; % deg
inc = 42; % deg
omega = 22; % deg (arg. of perigee)
nu = 0; % deg (true anomaly aka theta)

[r0, v0] = COEs_to_state([h, ecc, inc, RAAN, omega, nu], mu);
r0 = r0'; v0 = v0';

a=h^2/mu/(1 - ecc^2);
period = (2*pi*sqrt(a^3/mu));
n_sc = (2*pi)/period; % mean motion, rad/s

T_d0 = 10^-5 * [0; 0.5; 0]; % N*m
q_c = [0; 0; 0; 1]; % body-ECI

% LVLH-ECI
z0_lvlh = -r0 / norm(r0);
y0_lvlh = -(cross(r0,v0))/norm(cross(r0, v0));
x0_lvlh = cross(y0_lvlh, z0_lvlh);
C0_LVLH_ECI = [x0_lvlh'; y0_lvlh'; z0_lvlh'];
q0_LVLH_ECI = C2quat(C0_LVLH_ECI);

% Body-ECI
q0_body_ECI = [0.5; 0.5; 0.5; 0.5];
C0_body_ECI = quat2C(q0_body_ECI);
omega0_LVLH_ECI = C0_body_ECI*cross(r0, v0)/norm(r0)^2; % F_b, rad/s
omega0_body_ECI = 10^-3 * [0.5; -0.2; -3]; % F_b, rad/s
E0_body_ECI = euler_from_C(C0_body_ECI); % initial euler angles, rad

% Body-LVLH
C0_body_LVLH = C0_body_ECI * C0_LVLH_ECI';
q0_body_LVLH = C2quat(C0_body_LVLH);
omega0_body_LVLH = omega0_body_ECI - omega0_LVLH_ECI;
E0_body_LVLH = euler_from_C(C0_body_LVLH); % initial euler angles, rad

% Inertia
z = 3; % m
x = 1.5; y = x; % m
m_sc = 500; % kg

I_b = 1/12*m_sc*diag([y^2 + z^2; x^2 + z^2; x^2 + y^2]); % platform

% Control law
zeta = 0.7;
t_s = 70; % 2% settling time, s
omega_n = 4.4 / (t_s*zeta); % natural frequency, hz

CMG = boolean(1);
RW = boolean(0);

%% 3.1

theta = 57;

Ag0 = [sind(theta) 0 -sind(theta) 0; 0 sind(theta) 0 -sind(theta); cosd(theta) cosd(theta) cosd(theta) cosd(theta)];
At0 = [-cosd(theta) 0 cosd(theta) 0; 0 -cosd(theta) 0 cosd(theta); sind(theta) sind(theta) sind(theta) sind(theta)];

for i = 1:length(Ag0)
    As0(:,i) = cross(At0(:,i), Ag0(:,i));
end


disp("A_s = ")
disp(As0)
disp(" ")
disp("A_g = ")
disp(Ag0)
disp(" ")
disp("A_t = ")
disp(At0)

%% 3.2
m_CMG = 4.5; % kg
r_CMG = .2; % m
h_CMG = 5/100; % cm

I_xx = (1/12) * m_CMG * (h_CMG^2 + 3*r_CMG^2);
I_yy = I_xx;
I_zz = 0.5*m_CMG*r_CMG^2;

I_CMG1 = [I_xx 0 0; 0 I_yy 0; 0 0 I_zz]; % 1,1: I_s   2,2: I_t    3,3: I_g
I_CMG2 = I_CMG1; I_CMG3 = I_CMG1; I_CMG4 = I_CMG1;

I_ws = [I_CMG1(3,3) 0 0 0; 0 I_CMG2(3,3) 0 0; 0 0 I_CMG3(3,3) 0; 0 0 0 I_CMG4(3,3)];
I_ct = [I_CMG1(2,2) 0 0 0; 0 I_CMG2(2,2) 0 0; 0 0 I_CMG3(2,2) 0; 0 0 0 I_CMG4(2,2)];
I_cg = [I_CMG1(1,1) 0 0 0; 0 I_CMG2(1,1) 0 0; 0 0 I_CMG3(1,1) 0; 0 0 0 I_CMG4(1,1)];

I_cs = I_ws;

J_sc = I_b + As0*I_cs*As0' + At0*I_ct*At0' + Ag0*I_cg*Ag0';

disp("I_b = ")
disp(I_b)
disp(" ")
disp("J = ")
disp(J_sc)

%% 3.3
K_p = 2*omega_n^2 * J_sc;
K_d = 2*zeta*omega_n*J_sc;
disp('K_p = ')
disp(' ')
disp(K_p)
disp(" ")
disp('K_d = ')
disp(' ')
disp(K_d)


%% 3.4

Omega0 = [800; 800; 800; 800]; % rad/s
gamma0 = [0;0;0;0];

sim_time = 200;
simout = sim('HW3_sim.slx');

% q_b_ECI
figure
subplot(2,2,1)
plot(simout.tout, squeeze(simout.q_b_ECI))
legend('$\epsilon_1$', '$\epsilon_2$', '$\epsilon_3$', '$\eta$','fontsize', 8 ,'Interpreter', 'latex')
xlabel('time (sec)','fontsize', 8 ,'Interpreter', 'latex')
ylabel('Quaternion Component','fontsize', 8 ,'Interpreter', 'latex')
grid on
title('Quaternion from ${\mathcal{F}}_b$ to ${\mathcal{F}}_{ECI}$','fontsize', 12 ,'Interpreter', 'latex')
xlim([0 200])
set(gcf, 'Position', [100, 100, 800, 500]);

% E_b_ECI
subplot(2,2,2)
plot(simout.tout, squeeze(simout.E_b_ECI))
legend('$\phi$', '$\theta$', '$\psi$','fontsize', 8 ,'Interpreter', 'latex')
xlabel('time (sec)','fontsize', 8 ,'Interpreter', 'latex')
ylabel('Euler angle (rad)','fontsize', 8 ,'Interpreter', 'latex')
grid on
title('Euler Angles from ${\mathcal{F}}_b$ to ${\mathcal{F}}_{ECI}$','fontsize', 12 ,'Interpreter', 'latex')
xlim([0 200])
set(gcf, 'Position', [100, 100, 800, 500]);

% w_b_ECI
subplot(2,2,[3,4])
plot(simout.tout, squeeze(simout.w_b_ECI))
legend('$\omega_x$', '$\omega_y$', '$\omega_z$','fontsize', 8 ,'Interpreter', 'latex')
xlabel('time (sec)','fontsize', 8 ,'Interpreter', 'latex')
ylabel('Angular velocity (rad/s)','fontsize', 8 ,'Interpreter', 'latex')
grid on
title('$\omega_{b-ECI}$ in ${\mathcal{F}}_b$','fontsize', 12 ,'Interpreter', 'latex')
xlim([0 200])
set(gcf, 'Position', [100, 100, 800, 500]);

% Gimbal angles
figure
plot(simout.tout, squeeze(simout.gamma))
legend('$\gamma_1$', '$\gamma_2$', '$\gamma_3$', '$\gamma_4$', 'fontsize', 10 ,'Interpreter', 'latex')
xlabel('time (sec)','fontsize', 12 ,'Interpreter', 'latex')
ylabel('Gimbal angles (rads)','fontsize', 12 ,'Interpreter', 'latex')
grid on
title('Gimbal Angles','fontsize', 12 ,'Interpreter', 'latex')
xlim([0 200])


























%% Functions and such

function [r, v] = COEs_to_state(COEs, mu)

h = COEs(1);
e = COEs(2);
inc = COEs(3);
RAAN = COEs(4);
omega = COEs(5);
theta = COEs(6);

% h = sqrt(a*mu*(1-e^2));

rp = (h^2/mu) * (1/(1 + e*cosd(theta))) * (cosd(theta)*[1;0;0] + sind(theta)*[0;1;0]);
vp = (mu/h) * (-sind(theta)*[1;0;0] + (e + cosd(theta))*[0;1;0]);


Q = Qfunc(omega, inc, RAAN);
Q = Q';

r = (Q*rp)';
v = (Q*vp)';

end

function Q = Qfunc(omega, inc, RAAN)

R3_omega = [cosd(omega) sind(omega) 0; -sind(omega) cosd(omega) 0; 0 0 1];
R1_inc = [1 0 0; 0 cosd(inc) sind(inc); 0 -sind(inc) cosd(inc)];
R3_RAAN = [cosd(RAAN) sind(RAAN) 0; -sind(RAAN) cosd(RAAN) 0; 0 0 1];

Q = R3_omega * R1_inc * R3_RAAN;

end

function quat = C2quat(C)

eta = sqrt(trace(C) + 1)/2;

if eta == 0
    abs_eps_1 = sqrt(0.5*(C(1,1)+1));
    abs_eps_2 = sqrt(0.5*(C(2,2)+1));
    abs_eps_3 = sqrt(0.5*(C(3,3)+1));

    if abs_eps_1 > 0
        eps_1 = abs(abs_eps_1);
        eps_2 = sign(C(1,2)) * abs_eps_2;
        eps_3 = sign(C(1,3)) * abs_eps_3;
    elseif abs_eps_2 > 0
        eps_1 = sign(C(1,2)) * abs_eps_1;
        eps_2 = abs_eps_2;
        eps_3 = sign(C(2,3)) * abs_eps_3;
    else 
        eps_1 = sign(C(1,3)) * abs_eps_1;
        eps_2 = sign(C(2,3)) * abs_eps_2;
        eps_3 = abs_eps_3;
    end

else
    eps_1 = (C(2,3) - C(3,2))/(4*eta);
    eps_2 = (C(3,1) - C(1,3))/(4*eta);
    eps_3 = (C(1,2) - C(2,1))/(4*eta);
end

quat = [eps_1; eps_2; eps_3; eta];

end

function C = quat2C(quat)

C = ((2*(quat(4)^2) - 1)*eye(3)) + (2*quat(1:3)*...
    quat(1:3)') - (2*quat(4)*[0 -quat(3) ...
    quat(2); quat(3) 0 -quat(1); -quat(2) ...
    quat(1) 0]);

end

function angles = euler_from_C(C)

phi = atan2(C(2,3), C(3,3)); % rad
theta = -asin(C(1,3)); % rad
psi = atan2(C(1,2), C(1,1)); % rad

angles = [phi; theta; psi];

end
