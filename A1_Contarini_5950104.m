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

% Set proportional and derivative gains
kp1 = 5;
kp2 = 5;
kp3 = 5;
kd1 = 100;
kd2 = 100;
kd3 = 100;

kp = [kp1; kp2; kp3];
kd = [kd1; kd2; kd3];

% Integrate dynamics and kinematics equation
t_start = 0;
t_integration = 3*60*60;
dt = 0.1;
t_span = linspace(t_start, t_start + t_integration, t_integration/dt);

[t, y] = ode45(@(t, y)fun(y, J, kp, kd, n, Md, H), t_span, y0);

% Plot theta angles over time
figure(1)
plot(t, rad2deg(y(:,1)), LineWidth=1, Color='blue', LineStyle='-')
hold on
plot(t, rad2deg(y(:,2)), LineWidth=1, Color='red', LineStyle='-')
plot(t, rad2deg(y(:,3)), LineWidth=1, Color='green', LineStyle='-')
yline(0.1, Color='black', LineWidth=1)
yline(-0.1, Color='black', LineWidth=1)
legend('\theta_1', '\theta_2', '\theta_3', '\theta = \pm 0.1 deg', fontsize=15)
xlabel('Integration time  [s]', FontSize=15)
ylabel('Angle  [deg]', fontsize=15)
ax = gca(figure(1));
ax.FontSize = 15;
grid("on")
title('Performance of PD controller (ideal measurements)', FontSize=15)
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
plot(t2, rad2deg(y2(:,1)), LineWidth=1, Color='blue')
hold on
plot(t2, rad2deg(y2(:,2)), LineWidth=1, Color='red')
plot(t2, rad2deg(y2(:,3)), LineWidth=1, Color='green')
yline(0.1, Color='black', LineWidth=1)
yline(-0.1, Color='black', LineWidth=1)
legend('\theta_1', '\theta_2', '\theta_3', '\theta = \pm 0.1 deg', fontsize=15)
xlabel('Integration time  [s]', FontSize=15)
ylabel('Angle  [deg]', fontsize=15)
ax = gca(figure(2));
ax.FontSize = 15;
grid("on")
title('Performance of PD controller (with error measurements)', FontSize=15)
saveas(figure(2), 'task.4.pdf')
hold off

%% Task 7 - corrected


%% Task 7

dt = 10;
t2 = linspace(t_start, t_start + t_integration, t_integration/dt);

N = length(t2);

b_1_0 = deg2rad(5);
b_2_0 = deg2rad(5);
b_3_0 = deg2rad(5);

P_0 = diag([1 1 1 1 1 1]);

Phi_0 = eye(6);
Phi_0 = reshape(Phi_0, 36, 1);

P_k = P_0;

Q = diag([0.1, 0.1, 0.1, 0.1, 0.1, 0.1]);

y_k = [theta_1_0; theta_2_0; theta_3_0; omega_1_0; omega_2_0; omega_3_0; Hw_1_0; Hw_2_0; Hw_3_0];

R = diag((deg2rad([0.1 0.1 0.1 0.1 0.1 0.1])).^2);

theta_store = zeros(N-1, 3);
b_store = zeros(N-1, 3);

H_matrix = eye(6);

for k = 1:N-1
    
    % Store current state
    theta_store(k, :) = (y_k(1:3)).';
    
    % Generate attitude measurements errors
    attitude_measurement_error = normrnd(0, deg2rad(0.1), 1, 3);
    
    % Compute control moment of reaction wheels and control moment gyro
    theta_1 = y_k(1);
    theta_2 = y_k(2);
    theta_3 = y_k(3);
    omega_1 = y_k(4);
    omega_2 = y_k(5);
    omega_3 = y_k(6);
    theta_dot = compute_theta_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, n);
    Hw_dot_1 = kp(1) * theta_1 + kd(1) * theta_dot(1);
    Hw_dot_2 = kp(2) * theta_2 + kd(2) * theta_dot(2);
    Hw_dot_3 = kp(3) * theta_3 + kd(3) * theta_dot(3);
    
    % Integrate equations of motion and state transition matrix
    t_span = t2(k:k+1);
    [t3, y_k1_k] = ode45(@(t, y)fun3(y, J, n, Md, H, Hw_dot_1, Hw_dot_2, Hw_dot_3), t_span, y_k);
    
    % Add measurements errors and noise to produce the measurements
    y_k1_k = (y_k1_k(length(t3), :))';
    x_k1_k = y_k1_k(1:6);
    theta_1_measured = y_k1_k(1) + attitude_measurement_error(1);
    theta_2_measured = y_k1_k(2) + attitude_measurement_error(2);
    theta_3_measured = y_k1_k(3) + attitude_measurement_error(3);
    omega_measured = y_k1_k(4:6) + angular_velocity_measurement_errors;
    
    % Store measurements in one vector
    z_k1 = [theta_1_measured; theta_2_measured; theta_3_measured; omega_measured];

    % Compute predicted state
    [t4, predicted_state] = ode45(@(t, y)integrate_state_and_state_transion_matrix(y, J, kp, kd, n, Md, H), t_span, [z_k1; y_k(7); y_k(8); y_k(9); reshape(eye(6, 6), 36, 1)]);
    
    predicted_state = (predicted_state(length(t4), :))';

    % Reshape state transition matrix
    Phi_k1_k = reshape(predicted_state(10:45), 6, 6);
    
    % Propagate covariance matrix
    P_k1_k = Phi_k1_k * P_k * Phi_k1_k.' + Q;
    
    % Calculate Kalman gain matrix
    K_k1_k = P_k1_k * H_matrix.' /(H_matrix * P_k1_k * H_matrix.' + R);

    % Calculate predicted observations
    z_predicted = predicted_state(1:6);

    % Update state
    estimated_state = predicted_state(1:6) + K_k1_k * (z_k1 - z_predicted);
    
    % Update covariance matrix
    P_k1_k1 = (eye(6) - K_k1_k*H_matrix) * P_k1_k;

    % Prepare for next iteration
    y_k = [estimated_state; predicted_state(7:9)];
    P_k = P_k1_k1;
    
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
ax = gca(figure(3));
ax.FontSize = 15;
grid("on")
title('Performance of PD controller', FontSize=15)
saveas(figure(3), 'task.7.pdf')
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

theta_dot = compute_theta_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, n);

theta_dot_1 = theta_dot(1);
theta_dot_2 = theta_dot(2);
theta_dot_3 = theta_dot(3);

Hw_dot_1 = kp(1) * theta_1 + kd(1) * theta_dot_1;
Hw_dot_2 = kp(2) * theta_2 + kd(2) * theta_dot_2;
Hw_dot_3 = kp(3) * theta_3 + kd(3) * theta_dot_3;

omega_dot = compute_omega_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, J, Hw_dot_1, Hw_dot_2, Hw_dot_3, Hw_1, H, Hw_3, Md, n);

omega_dot_1 = omega_dot(1);
omega_dot_2 = omega_dot(2);
omega_dot_3 = omega_dot(3);

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

% Function for the prediction step of the state and of the state transition
% matrix
function dy = fun3(y, J, n, Md, H, Hw_dot_1, Hw_dot_2, Hw_dot_3)

theta_1 = y(1);
theta_2 = y(2);
theta_3 = y(3);
omega_1 = y(4);
omega_2 = y(5);
omega_3 = y(6);
Hw_1 = y(7);
Hw_2 = y(8);
Hw_3 = y(9);

theta_dot = compute_theta_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, n);

theta_dot_1 = theta_dot(1);
theta_dot_2 = theta_dot(2);
theta_dot_3 = theta_dot(3);

omega_dot = compute_omega_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, J, Hw_dot_1, Hw_dot_2, Hw_dot_3, Hw_1, H, Hw_3, Md, n);

omega_dot_1 = omega_dot(1);
omega_dot_2 = omega_dot(2);
omega_dot_3 = omega_dot(3);

dy = [theta_dot_1; theta_dot_2; theta_dot_3; omega_dot_1; omega_dot_2; omega_dot_3; Hw_dot_1; Hw_dot_2; Hw_dot_3];

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

function F = compute_F_matrix(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, n, J, kp, H, Hw_1, Hw_3)

F = zeros(6, 6);
F(1, 1) = omega_2 * cos(theta_1) * tan(theta_2) - omega_3 * sin(theta_1)*tan(theta_2);
F(1, 2) = omega_2 * sin(theta_1) / (cos(theta_2))^2 + omega_3 * cos(theta_1) / (cos(theta_2))^2 + n * sin(theta_2) * sin(theta_3) / (cos(theta_2))^2;
F(1, 3) = n * cos(theta_3) / cos(theta_3);
F(1, 4) = 1;
F(1, 5) = sin(theta_1) * tan(theta_2);
F(1, 6) = cos(theta_1) * tan(theta_2);

F(2, 1) = - sin(theta_1) * omega_2 - cos(theta_1) * omega_3;
F(2, 3) = - n * sin(theta_3);
F(2, 5) = cos(theta_1);
F(2, 6) = -sin(theta_1);

F(3, 1) = omega_2 * cos(theta_1) / cos(theta_2) - omega_3 * sin(theta_1) / cos(theta_2);
F(3, 2) = omega_2 * sin(theta_1) * sin(theta_2) / (cos(theta_2))^2 - omega_3 * sin(theta_1) * sin(theta_2) / (cos(theta_2))^2;
F(3, 3) = n * sin(theta_2) * cos(theta_3) / cos(theta_2);
F(3, 5) = sin(theta_1)/cos(theta_2);
F(3, 6) = cos(theta_1)/cos(theta_2);

F(4, 1) = (1/J(1, 1)) * ( -kp(1) -3 * n^2 * (J(2, 2) - J(3, 3)) * ((cos(theta_1))^2 * (cos(theta_2))^2 - (sin(theta_1))^2 * (cos(theta_2))^2 ));
F(4, 2) = (1/J(1, 1)) * 6 * n^2 * (J(2, 2) - J(3, 3)) * sin(theta_1) * cos(theta_1) * cos(theta_2) * sin(theta_2);
F(4, 5) = (1/J(1, 1)) * ((J(2, 2) - J(3,3))*omega_3 - Hw_3);
F(4, 6) = (1/J(1, 1)) * ((J(2, 2) - J(3, 3))*omega_2 - H);

F(5, 1) = (1/J(2, 2)) * (-3 * n^2 * (J(3, 3) - J(1, 1)) * sin(theta_2) * sin(theta_1) * cos(theta_2));
F(5, 2) = (1/J(2, 2)) * (-3 * n^2 * (J(3, 3) - J(1, 1)) * (-(cos(theta_2))^2 * cos(theta_1) + (sin(theta_2))^2 * cos(theta_1)) - kp(2));
F(5, 4) = (J(3, 3) - J(1, 1)) * omega_3 / J(2, 2);
F(5, 6) = (J(3, 3) - J(1, 1)) * omega_1 / J(2, 2);

F(6, 1) = (1/J(3, 3)) * 3 * n^2 * (J(1, 1) - J(2, 2)) * sin(theta_2) * cos(theta_1) * cos(theta_2);
F(6, 2) = (1/J(3, 3)) * (-3 * n^2 * (J(1, 1) - J(2, 2)) * (-cos(theta_2)^2 * cos(theta_1) + sin(theta_2)^2 * cos(theta_1)));
F(6, 3) = - kp(3) / J(3, 3);
F(6, 4) = (H + omega_2*(J(1, 1) - J(2, 2))) / J(3, 3);
F(6, 5) = (1/J(3, 3))*(Hw_1 + (J(1, 1) - J(2, 2))*omega_1);

end

function dy = integrate_state_and_state_transion_matrix(y, J, kp, kd, n, Md, H)

theta_1 = y(1);
theta_2 = y(2);
theta_3 = y(3);
omega_1 = y(4);
omega_2 = y(5);
omega_3 = y(6);
Hw_1 = y(7);
Hw_2 = y(8);
Hw_3 = y(9);
Phi = reshape(y(10:45), 6, 6);

theta_dot = compute_theta_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, n);

theta_dot_1 = theta_dot(1);
theta_dot_2 = theta_dot(2);
theta_dot_3 = theta_dot(3);

Hw_dot_1 = kp(1) * theta_1 + kd(1) * theta_dot_1;
Hw_dot_2 = kp(2) * theta_2 + kd(2) * theta_dot_2;
Hw_dot_3 = kp(3) * theta_3 + kd(3) * theta_dot_3;

omega_dot = compute_omega_dot_vector(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, J, Hw_dot_1, Hw_dot_2, Hw_dot_3, Hw_1, H, Hw_3, Md, n);

omega_dot_1 = omega_dot(1);
omega_dot_2 = omega_dot(2);
omega_dot_3 = omega_dot(3);

F = compute_F_matrix(theta_1, theta_2, theta_3, omega_1, omega_2, omega_3, n, J, kp, H, Hw_1, Hw_3);
Phi_dot = F*Phi;
Phi_dot = reshape(Phi_dot, 36, 1);

dy = [theta_dot_1; theta_dot_2; theta_dot_3; omega_dot_1; omega_dot_2; omega_dot_3; Hw_dot_1; Hw_dot_2; Hw_dot_3; Phi_dot];

end
