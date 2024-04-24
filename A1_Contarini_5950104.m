% Script for Assignment 1 of Spacecraft Attitude Dynamics and Control -
% Mattia Contarini (5950104)

clear
clc
close all

%% Constants and given parameters

gravitational_parameter = 3.986004418 * 1e14;  % Earth gravitational 
% parameter  [m^3/s^3]

Re = 6378.137 * 1e3;  % Earth equatorial radius [m]
h = 700 * 1e3;  % Circular orbit altitude [m]

r = Re + h;

n = sqrt(gravitational_parameter / r^3);

J = [2500 0 0;
    0 2300 0;
    0 0 3000];

Md = [0.001; 0.001; 0.001];

theta_ss = deg2rad([0.1; 0.1; 0.1]); % Required steady-state error

damping_ratio = [0.707; 0.707; 0.707];

H = 38.2; % Angular momentum of momentum bias wheel [N*m*s]. Retrieved from
% SMAD, Table 11-12

%% Task 3
% Simulate the performance of the controller with perfect attitude and
% angular velocity measurements without sensor noise and bias

% Compute proportional and derivative control parameters for momentum bias
% wheel on first, second and thid axis
kp2 = Md(2) / theta_ss(2);
kd2 = 2 * damping_ratio(2) * sqrt(kp2 * J(2, 2));
kp1 = (1 - n*H*theta_ss(1)/Md(1)) / (theta_ss(1) / Md(1));
kp3 = (1 - n*H*theta_ss(3)/Md(3)) / (theta_ss(3) / Md(3));
u = sqrt((kp1 * kp3 + (n*H)^2 + n*H*(kp1 + kp3)) / (J(1, 1) * J(3, 3)));
v = H^2 + J(3, 3) * (kp1 + n * H) + J(1, 1) * (kp3 + n * H);
w = (-kp3 - n * H + u * J(3, 3)) / (kp1 + n * H - u * J(1, 1));
kd1 = sqrt((v - u * J(1, 1) * J(3, 3) * (4 * damping_ratio(2)^2 - 2)) / ...
    ((J(1, 1) / (4 * damping_ratio(2)^2 * J(1, 1))) + ...
    (w / (2 * damping_ratio(2)^2)) + ...
    ((w^2 * J(1, 1)) / (4 * damping_ratio(2)^2 * J(3, 3))) - w));
kd3 = w * kd1;

% Set initial elements of state vector y
theta_1_0 = deg2rad(10);
theta_2_0 = deg2rad(10);
theta_3_0 = deg2rad(10);
omega_1_0 = - n * theta_3_0;
omega_2_0 = -n;
omega_3_0 = n * theta_1_0;
Hw_1_0 = 0;
Hw_2_0 = 0;
Hw_3_0 = 0;

% Initialize state vector
y0 = [theta_1_0 theta_2_0 theta_3_0 omega_1_0 omega_2_0 omega_3_0 Hw_1_0 Hw_2_0 Hw_3_0];

kp = [kp1; kp2; kp3];
kd = [kd1; kd2; kd3];

% Integrate dynamics and kinematics equation
t_start = 0;
t_integration = 20*60;
dt = 0.1;
t_span = linspace(t_start, t_start + t_integration, t_integration/dt);

[t, y] = ode45(@(t, y)fun(y, J, kp, kd, n, Md, H), t_span, y0);

% Plot theta angles over time
figure(1)
plot(t, rad2deg(y(:,1)), LineWidth=2, Color='blue')
hold on
plot(t, rad2deg(y(:,2)), LineWidth=2, Color='red')
plot(t, rad2deg(y(:,3)), LineWidth=2, Color='green')
yline(0.1, Color='black', LineWidth=2)
legend('\theta_1', '\theta_2', '\theta_3', '\theta = 0.1 deg', fontsize=15)
xlabel('Integration time  [s]', FontSize=15)
ylabel('Angle  [deg]', fontsize=15)
ax = gca(figure(1));
ax.FontSize = 15;
grid("on")
title('Performance of PD controller', FontSize=15)
saveas(figure(1), 'task.3.pdf')
hold off

%% Task 4

% Set initial elements of state vector y
attitude_measurement_error = normrnd(0, deg2rad(0.1));
theta_1_0 = deg2rad(10) - attitude_measurement_error;
theta_2_0 = deg2rad(10) - attitude_measurement_error;
theta_3_0 = deg2rad(10) - attitude_measurement_error;
angular_velocity_measurement_errors = deg2rad([0.1 -0.1 0.15]);
omega_1_0 = - n * theta_3_0 - angular_velocity_measurement_errors(1);
omega_2_0 = -n - angular_velocity_measurement_errors(2);
omega_3_0 = n * theta_1_0 - angular_velocity_measurement_errors(3);
Hw_1_0 = 0;
Hw_2_0 = 0;
Hw_3_0 = 0;

% Initialize state vector
y0 = [theta_1_0 theta_2_0 theta_3_0 omega_1_0 omega_2_0 omega_3_0 Hw_1_0 Hw_2_0 Hw_3_0];

[t2, y2] = ode45(@(t, y)fun(y, J, kp, kd, n, Md, H), t_span, y0);

% Plot theta angles over time
figure(2)
plot(t2, rad2deg(y2(:,1)), LineWidth=2, Color='blue')
hold on
plot(t2, rad2deg(y2(:,2)), LineWidth=2, Color='red')
plot(t2, rad2deg(y2(:,3)), LineWidth=2, Color='green')
yline(0.1, Color='black', LineWidth=2)
legend('\theta_1', '\theta_2', '\theta_3', '\theta = 0.1 deg', fontsize=15)
xlabel('Integration time  [s]', FontSize=15)
ylabel('Angle  [deg]', fontsize=15)
ax = gca(figure(2));
ax.FontSize = 15;
grid("on")
title('Performance of PD controller', FontSize=15)
saveas(figure(2), 'task.4.pdf')
hold off


%% Functions

function dy = fun(y, J, kp, kd, n, Md, H)
theta_1 = y(1);
theta_2 = y(2);
theta_3 = y(3);
omega_1 = y(4);
omega_2 = y(5);
omega_3 = y(6);
Hw_1 = y(7);
Hw_2 = y(8);
Hw_3 = y(9);

theta_dot_1 = omega_1 + n * theta_3;
theta_dot_2 = omega_2 + n;
theta_dot_3 = omega_3 - n * theta_1;

Hw_dot_1 = kp(1) * theta_1 + kd(1) * theta_dot_1;
Hw_dot_2 = kp(2) * theta_2 + kd(2) * theta_dot_2;
Hw_dot_3 = kp(3) * theta_3 + kd(3) * theta_dot_3;

omega_dot_1 = (1/J(1, 1)) * (Md(1) - 3 * n^2 * (J(2, 2) - J(3, 3)) * theta_1 - (J(2, 2) - J(3, 3)) * n * omega_3 - Hw_dot_1 - omega_3 * H - n*Hw_3);
omega_dot_2 = (1/J(2, 2)) * (3 * n^2 *(J(3, 3) - J(1, 1)) * theta_2 + Md(2) - Hw_dot_2);
omega_dot_3 = (1/J(3, 3)) * (Md(3) + omega_1 * H + n * Hw_1 - Hw_dot_3 - (J(1, 1) - J(2, 2))*n*omega_1);

dy = [theta_dot_1; theta_dot_2; theta_dot_3; omega_dot_1; omega_dot_2; omega_dot_3; Hw_dot_1; Hw_dot_2; Hw_dot_3];

end
