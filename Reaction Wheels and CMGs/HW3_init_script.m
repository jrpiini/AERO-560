%% AERO 560 HW 3 - Joseph Piini

clear 
close all
clc

addpath('C:\Users\josep\OneDrive - Cal Poly\Desktop\AERO 560\Functions')

%% Problem 1

mu = 398600.4418;
alt = 700;
r_earth = 6378;

speed = sqrt(mu/(alt+r_earth)); % km/s
period = 2*pi*sqrt((alt+r_earth)^3 / mu);
% [~, rho] = atmosnrlmsise00(alt*1000, 35, -120, 2025, 38, 0); % kg/m^3
rho = 3.91e-14; % MSISE-90, mean solar activity
Cd = 2.2;
A_wetted = 3;

F_drag = 0.5*rho*(speed*1000)^2*Cd*A_wetted;
torque_drag = F_drag*0.25;
momentum = torque_drag*period;

disp("Assumptions/eqns: C_d = 2.2, rho = 3.91e-14 kg/m^3 (MSISE-90, mean solar activity), F_drag = 0.5*rho*v^2*Cd*A, torque = F_drag*0.25")
disp(" ")
disp("Momentum accumulated over orbit = " + momentum + " Nms")
disp(" ")
disp("The Rocket Lab RW-0.03 would be sufficient as it can store 30 mNms of momentum, roughly 3 times the amount accumulated over an orbit. The spec sheet is at the end of the pdf.")

%% Problem 2

CMG = boolean(0);
RW = boolean(1);

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

% Inertia
z = 3; % m
x = 1.5; y = x; % m
m_sc = 500; % kg

I_b = 1/12*m_sc*diag([y^2 + z^2; x^2 + z^2; x^2 + y^2]); % platform

% Control law
zeta = 0.7;
t_s = 70; % 2% settling time, s
omega_n = 4.4 / (t_s*zeta); % natural frequency, hz

%% 2.1

% Wheel pyramid
theta = 57; % wheel inclination

As0 = [sind(theta) 0 -sind(theta) 0; 0 sind(theta) 0 -sind(theta); cosd(theta) cosd(theta) cosd(theta) cosd(theta)];
At0 = [-cosd(theta) 0 cosd(theta) 0; 0 -cosd(theta) 0 cosd(theta); sind(theta) sind(theta) sind(theta) sind(theta)];

for i = 1:length(As0)
    Ag0(:,i) = cross(At0(:,i), As0(:,i));
end

disp("A_s = ")
disp(As0)
disp(" ")
disp("A_g = ")
disp(Ag0)
disp(" ")
disp("A_t = ")
disp(At0)

%% 2.2
m_wheel = 2; % kg
r_wheel = .08; % m
h_wheel = .02; % m

I_xx = (1/12) * m_wheel * (h_wheel^2 + 3*r_wheel^2);
I_yy = I_xx;
I_zz = 0.5*m_wheel*r_wheel^2;

I_wheel1 = [I_xx 0 0; 0 I_yy 0; 0 0 I_zz]; % 1,1: I_s   2,2: I_t    3,3: I_g
I_wheel2 = I_wheel1; I_wheel3 = I_wheel1; I_wheel4 = I_wheel1;

I_cg = [I_wheel1(1,1) 0 0 0; 0 I_wheel2(1,1) 0 0; 0 0 I_wheel3(1,1) 0; 0 0 0 I_wheel4(1,1)];
I_ct = [I_wheel1(2,2) 0 0 0; 0 I_wheel2(2,2) 0 0; 0 0 I_wheel3(2,2) 0; 0 0 0 I_wheel4(2,2)];
I_ws = [I_wheel1(3,3) 0 0 0; 0 I_wheel2(3,3) 0 0; 0 0 I_wheel3(3,3) 0; 0 0 0 I_wheel4(3,3)];

I_cs = I_ws;

J_sc = I_b + As0*I_cs*As0' + At0*I_ct*At0' + Ag0*I_cg*Ag0';

disp("I_b = ")
disp(I_b)
disp(" ")
disp("J = ")
disp(J_sc)

%% 2.3
K_p = 2*omega_n^2 * J_sc;
K_d = 2*zeta*omega_n*J_sc;
disp('K_p = ')
disp(' ')
disp(K_p)
disp(" ")
disp('K_d = ')
disp(' ')
disp(K_d)


%% 2.4

T_d0 = 10^-5 * [0; 0.5; 0]; % N*m
q_c = [0; 0; 0; 1];

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
% C0_body_LVLH = eye(3); % initially aligned
q0_body_LVLH = C2quat(C0_body_LVLH);
omega0_body_LVLH = omega0_body_ECI - omega0_LVLH_ECI;
E0_body_LVLH = euler_from_C(C0_body_LVLH); % initial euler angles, rad

% Reaction wheels
Omega0 = [0;0;0;0]; % initial speeds
gamma0 = [0;0;0;0]; % no CMGs yet

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

%% 2.5

% RW Speeds
figure
plot(simout.tout, squeeze(simout.RW_speeds))
legend('$\Omega_1$', '$\Omega_2$', '$\Omega_3$', '$\Omega_4$', 'fontsize', 10 ,'Interpreter', 'latex')
xlabel('time (sec)','fontsize', 12 ,'Interpreter', 'latex')
ylabel('Angular velocity (rad/s)','fontsize', 12 ,'Interpreter', 'latex')
grid on
title('$\Omega$ in wheel space spin directions','fontsize', 12 ,'Interpreter', 'latex')
xlim([0 200])


%% 2.6

disp("We can see from the plots that the reaction wheels initially speed up to bring the spacecraft to the commanded quaternion, before settling to steady-state values after roughly 70 seconds (settling time). After this, each wheel maintains roughly the same speed ensuring the spacecraft retains the correct attitude with near-zero angular body rates.")

%% 2.7

disp("Matlab init script and simulink model are attached to submission")

%% Problem 3

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
clear simout
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

%% 3.5

% Gimbal angles
figure
plot(simout.tout, squeeze(simout.gamma))
legend('$\gamma_1$', '$\gamma_2$', '$\gamma_3$', '$\gamma_4$', 'fontsize', 10 ,'Interpreter', 'latex')
xlabel('time (sec)','fontsize', 12 ,'Interpreter', 'latex')
ylabel('Gimbal angles (rads)','fontsize', 12 ,'Interpreter', 'latex')
grid on
title('Gimbal Angles','fontsize', 12 ,'Interpreter', 'latex')
xlim([0 200])

%% 3.6

disp("We can see from the plots that the CMGs have nearly the same effect as the reaction wheels in that they are able to quickly and effectively lead the spacecraft from its initial state to the commanded attitude. There is a quick initial change in gimbal angles as the CMGs move, but shortly after 70 seconds, the gimbal angles remain pretty much constant and the spacecraft remains in the same orientation with near-zero body rates.")

%% 3.7

disp("Matlab init script and simulink model are attached to submission")