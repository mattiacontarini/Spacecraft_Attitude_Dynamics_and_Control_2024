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

angular_velocity_measurement_errors = deg2rad([0.1; -0.1; 0.15]);

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

[t2, y2] = ode45(@(t, y)fun2(y, J, kp, kd, n, Md, H), t_span, y0);

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


%% Task 7

N = length(t2);

b_1_0 = deg2rad(5);
b_2_0 = deg2rad(5);
b_3_0 = deg2rad(5);

P_0 = diag([100 100 100 100 100 100]);

Phi_0 = eye(6);
Phi_0 = reshape(Phi_0, 36, 1);

F = [0 0 n -1 0 0;
    0 0 0 0 -1 0;
    -n 0 0 0 0 -1;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];

H_matrix = eye(6);

P_k = P_0;

Q = zeros(6);

y_k = [theta_1_0; theta_2_0; theta_3_0; omega_1_0; omega_2_0; omega_3_0; Hw_1_0; Hw_2_0; Hw_3_0; b_1_0; b_2_0; b_3_0; Phi_0];

R = diag((deg2rad([0.1 0.1 0.1 0 0 0])).^2);

theta_store = zeros(N-1, 3);
b_store = zeros(N-1, 3);

for k = 1:N-1

    theta_store(k, :) = (y_k(1:3))';
    b_store(k, :) = (y_k(10:12))';

    attitude_measurement_error = normrnd(0, deg2rad(0.1), 1, 3);
    
    t_span = t2(k:k+1);
    [t3, y_k1_k] = ode45(@(t, y)fun3(y, J, kp, kd, n, Md, H, attitude_measurement_error, F), t_span, y_k);
    
    y_k1_k = (y_k1_k(length(t3), :))';

    x_k1_k = [y_k1_k(1:3); y_k1_k(10:12)];

    theta_1_measured = y_k1_k(1) + attitude_measurement_error(1);
    theta_2_measured = y_k1_k(2) + attitude_measurement_error(2);
    theta_3_measured = y_k1_k(3) + attitude_measurement_error(3);
    omega_measured = y_k1_k(4:6) + angular_velocity_measurement_errors;
   
    z_k1 = [theta_1_measured; theta_2_measured; theta_3_measured; omega_measured];

    Phi_k1_k = reshape(y_k1_k(13:48), 6, 6);

    P_k1_k = Phi_k1_k * P_k * Phi_k1_k' + Q;

    K_k1_k = P_k1_k * H_matrix' * inv(H_matrix * P_k1_k * H_matrix' + R);
    
    z_estimated = [y_k1_k(1:3) + attitude_measurement_error'; y_k1_k(4:6) + y_k1_k(10:12)];

    x_k1_k1 = x_k1_k + K_k1_k * (z_k1 - z_estimated);

    P_k1_k1 = (eye(6) - K_k1_k*H_matrix) * P_k1_k;
    
    y_k1_k(1:3) = x_k1_k1(1:3);
    y_k1_k(10:12) = x_k1_k1(4:6);
    
    y_k = y_k1_k;
    P_k1_k = P_k1_k1;
    


end

% Plot theta angles over time
figure(3)
dim = length(t2);
plot(t2(1:dim-1), rad2deg(theta_store(:, 1)), LineWidth=2, Color='blue')
hold on
plot(t2(1:dim-1), rad2deg(theta_store(:, 2)), LineWidth=2, Color='red')
plot(t2(1:dim-1), rad2deg(theta_store(:, 3)), LineWidth=2, Color='green')
yline(0.1, Color='black', LineWidth=2)
legend('\theta_1', '\theta_2', '\theta_3', '\theta = 0.1 deg', fontsize=15)
xlabel('Integration time  [s]', FontSize=15)
ylabel('Angle  [deg]', fontsize=15)
ax = gca(figure(2));
ax.FontSize = 15;
grid("on")
title('Performance of PD controller', FontSize=15)
saveas(figure(3), 'task.7.pdf')
hold off

figure(4)
dim = length(t2);
plot(t2(1:dim-1), rad2deg(b_store(:, 1)), LineWidth=2, Color='blue')
hold on
plot(t2(1:dim-1), rad2deg(b_store(:, 2)), LineWidth=2, Color='red')
plot(t2(1:dim-1), rad2deg(b_store(:, 3)), LineWidth=2, Color='green')
yline(0.1, Color='blue', LineStyle='--' , LineWidth=2)
yline(-0.1, Color='red', LineStyle='--', LineWidth=2)
yline(0.15, Color='green', LineStyle='--', LineWidth=2)
legend('b_1', 'b_2', 'b_3', 'b = 0.1 deg/s', 'b = -0.1 deg/s', 'b = 0.15 deg/s', fontsize=15)
xlabel('Integration time  [s]', FontSize=15)
ylabel('Angle  [deg]', fontsize=15)
ax = gca(figure(2));
ax.FontSize = 15;
grid("on")
title('Performance of PD controller', FontSize=15)
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


function dy = fun2(y, J, kp, kd, n, Md, H)
attitude_measurement_error = normrnd(0, deg2rad(0.1), 1, 3);

theta_1 = y(1);
theta_2 = y(2);
theta_3 = y(3);
omega_1 = y(4);
omega_2 = y(5);
omega_3 = y(6);
Hw_1 = y(7);
Hw_2 = y(8);
Hw_3 = y(9);

angular_velocity_measurement_errors = deg2rad([0.1, -0.1, 0.15]);

omega_1_measured = omega_1 + angular_velocity_measurement_errors(1);
omega_2_measured = omega_2 + angular_velocity_measurement_errors(2);
omega_3_measured = omega_3 + angular_velocity_measurement_errors(3);

theta_1_measured = theta_1 + attitude_measurement_error(1);
theta_2_measured = theta_2 + attitude_measurement_error(2);
theta_3_measured = theta_3 + attitude_measurement_error(3);

theta_dot = compute_theta_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, n);

% So here I use both the measured omega and measured theta to determine the measured theta_dot to give to the controller 
theta_dot_from_measurements = compute_theta_dot_vector(theta_1_measured, theta_2_measured, theta_3_measured, omega_1_measured, omega_2_measured, omega_3_measured, n);

%theta_dot_1 = omega_1 + angular_velocity_measurement_errors(1) + n * theta_3;
%theta_dot_2 = omega_2 + angular_velocity_measurement_errors(2) + n;
%theta_dot_3 = omega_3 + angular_velocity_measurement_errors(3) - n * theta_1;

theta_dot_1 = theta_dot(1);
theta_dot_2 = theta_dot(2);
theta_dot_3 = theta_dot(3);

Hw_dot_1 = kp(1) * theta_1_measured + kd(1) * theta_dot_from_measurements(1);
Hw_dot_2 = kp(2) * theta_2_measured + kd(2) * theta_dot_from_measurements(2);
Hw_dot_3 = kp(3) * theta_3_measured + kd(3) * theta_dot_from_measurements(3);

omega_dot = compute_omega_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, J, Hw_dot_1, Hw_dot_2, Hw_dot_3, Hw_1, H, Hw_3, Md, n);

omega_dot_1 = omega_dot(1);
omega_dot_2 = omega_dot(2);
omega_dot_3 = omega_dot(3);

dy = [theta_dot_1; theta_dot_2; theta_dot_3; omega_dot_1; omega_dot_2; omega_dot_3; Hw_dot_1; Hw_dot_2; Hw_dot_3];

end

function dy = fun3(y, J, kp, kd, n, Md, H, attitude_measurement_error, F)

theta_1 = y(1);
theta_2 = y(2);
theta_3 = y(3);
omega_1 = y(4);
omega_2 = y(5);
omega_3 = y(6);
Hw_1 = y(7);
Hw_2 = y(8);
Hw_3 = y(9);
b_1 = y(10);
b_2 = y(11);
b_3 = y(12);
Phi = reshape(y(13:48), 6, 6);

angular_velocity_measurement_errors = deg2rad([0.1, -0.1, 0.15]);

omega_measured_1 = omega_1 + angular_velocity_measurement_errors(1);
omega_measured_2 = omega_2 + angular_velocity_measurement_errors(2);
omega_measured_3 = omega_3 + angular_velocity_measurement_errors(3);

theta_dot_1 = omega_measured_1 - b_1 + n * theta_3;
theta_dot_2 = omega_measured_2 - b_2 + n;
theta_dot_3 = omega_measured_3 - b_3 - n * theta_1;

theta_1_measured = theta_1 + attitude_measurement_error(1);
theta_2_measured = theta_2 + attitude_measurement_error(2);
theta_3_measured = theta_3 + attitude_measurement_error(3);

Hw_dot_1 = kp(1) * theta_1_measured + kd(1) * theta_dot_1;
Hw_dot_2 = kp(2) * theta_2_measured + kd(2) * theta_dot_2;
Hw_dot_3 = kp(3) * theta_3_measured + kd(3) * theta_dot_3;

omega_dot_1 = (1/J(1, 1)) * (Md(1) - 3 * n^2 * (J(2, 2) - J(3, 3)) * theta_1 - (J(2, 2) - J(3, 3)) * n * omega_3 - Hw_dot_1 - omega_3 * H - n*Hw_3);
omega_dot_2 = (1/J(2, 2)) * (3 * n^2 *(J(3, 3) - J(1, 1)) * theta_2 + Md(2) - Hw_dot_2);
omega_dot_3 = (1/J(3, 3)) * (Md(3) + omega_1 * H + n * Hw_1 - Hw_dot_3 - (J(1, 1) - J(2, 2))*n*omega_1);

b_dot_1 = 0;
b_dot_2 = 0;
b_dot_3 = 0;

Phi_dot = F*Phi;
Phi_dot = reshape(Phi_dot, 36, 1);

dy = [theta_dot_1; theta_dot_2; theta_dot_3; omega_dot_1; omega_dot_2; omega_dot_3; Hw_dot_1; Hw_dot_2; Hw_dot_3; b_dot_1; b_dot_2; b_dot_3; Phi_dot];

end

function theta_dot = compute_theta_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, n)

v = [sin(theta_3); cos(theta_2)*cos(theta_3); sin(theta_2)*sin(theta_3)];

omega = [omega_1; omega_2; omega_3];

matrix = [cos(theta_2) sin(theta_1)*sin(theta_2) cos(theta_1)*sin(theta_2); 0 cos(theta_1)*cos(theta_2) -sin(theta_1)*cos(theta_2); 0 sin(theta_1) cos(theta_1)];

theta_dot = (1/cos(theta_2)) * matrix * omega + (n/cos(theta_2)) * v;

end

function omega_dot = compute_omega_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, J, Hw_dot_1, Hw_dot_2, Hw_dot_3, Hw_1, H, Hw_3, Md, n)

C = [cos(theta_2)*cos(theta_3) cos(theta_2)*sin(theta_3) -sin(theta_2);
    sin(theta_1)*sin(theta_2)*cos(theta_3)-cos(theta_1)*sin(theta_3) sin(theta_1)*sin(theta_2)*sin(theta_3)+cos(theta_1)*cos(theta_3) sin(theta_1)*cos(theta_2);
    cos(theta_1)*sin(theta_2)*cos(theta_3)+sin(theta_1)*sin(theta_3) cos(theta_1)*sin(theta_2)*sin(theta_3)-sin(theta_1)*cos(theta_3) cos(theta_1)*cos(theta_2)];

omega_dot_1 = (1/J(1, 1)) * (Md(1) - 3 * n^2 * (J(2, 2) - J(3, 3)) * C(2, 3) * C(3, 3) + (J(2, 2) - J(3, 3)) * omega_2 * omega_3 - Hw_dot_1 - omega_3 * H - omega_2*Hw_3);
omega_dot_2 = (1/J(2, 2)) * (-3 * n^2 * (J(3, 3) - J(1, 1)) * C(3, 3) * C(1, 3) + Md(2) - Hw_dot_2 + (J(3, 3) - J(1, 1)) * omega_3 * omega_1);
omega_dot_3 = (1/J(3, 3)) * (-3 * n^2 * (J(1, 1) - J(2, 2)) * C(1, 3) * C(2, 3) + Md(3) + omega_1 * H + omega_2 * Hw_1 - Hw_dot_3 + (J(1, 1) - J(2, 2))*omega_2*omega_1);

omega_dot = [omega_dot_1; omega_dot_2; omega_dot_3];
end

