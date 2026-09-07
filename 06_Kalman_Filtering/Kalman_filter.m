clc; clear; close all;

% read data
m = readtable('Measeurements.xlsx', 'Sheet','Measurements');
m = m(1:end,1:4);
m = renamevars(m,["Var1","Var2","Var3","Var4"],["Time_s_","e_m_","n_m_","speed_m_s_"]);
tv = readtable('Measeurements.xlsx', 'Sheet','True values');

% initial state
ve_0 = 3.53; % m/s
vn_0 = 0.86;

% given values
PSD = 0.01;  % m^2 s^-3
%PSD = 0.1;  % m^2 s^-3
%PSD = 0.05;  % m^2 s^-3
sigma_e = 3; % m
sigma_n = 3; % m
sigma_v = 0.5; % m/s
% sigma_e = 3*0.5; sigma_n = 3*0.5; sigma_v = 0.5*0.5;
% sigma_e = 3*0.01; sigma_n = 3*0.01; sigma_v = 0.5*0.01;
% sigma_e = 3*2; sigma_n = 3*2; sigma_v = 0.5*2;
sigma_v0_e = 3; % m/s
sigma_v0_n = 3; % m/s
sigma_0_e = 10; % m
sigma_0_n = 10; % m
%sigma_v0_e = 3*0.5; sigma_v0_n = 3*0.5; sigma_0_e = 10*0.5; sigma_0_n = 10*0.5; 
%sigma_v0_e = 3*2; sigma_v0_n = 3*2; sigma_0_e = 10*2; sigma_0_n = 10*2;



delta_t = 2; 

N = length(m.Time_s_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearised position-velocity model
% Step 1: Initialisation
F = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]; % system dynamic matrix
G = [0 0; 0 0; 1 0; 0 1]; % distribution matrix

% transition matrix
Tk = eye(4) + delta_t*F;
% covariance matrix of observations
R = diag([sigma_e^2,sigma_n^2,sigma_v^2]);
% covariance matrix of initial state vector
Q0 = diag([sigma_0_e^2,sigma_0_n^2,sigma_v0_e^2,sigma_v0_n^2]);
% process noise covariance matrix
Q = [PSD 0; 0 PSD];
Q_G = G*Q*G';
Qk = Q_G*delta_t + (F*Q_G + Q_G*F')*delta_t^2/2 + F*Q_G*F'*delta_t^3/3;

% initial state vector
x0 = [m.e_m_(1),m.n_m_(1),ve_0,vn_0]'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman filtering

x_pred = zeros(4,N-1); 
x_update = zeros (4,N); 
Q_pred = zeros(4,4,N-1); 
Q_update = zeros(4,4,N); 

x_update(:,1) = x0;
Q_update(:,:,1) = Q0;

x = x0;
Qx = Q0;

for k = 2: N
    % Step 2: time propagation
    x_pred(:,k) = Tk * x;
    Q_pred(:,:,k) = Tk * Qx * Tk' + Qk;

    % Step 3: Gain calculation
    v_ = sqrt(x_pred(3,k)^2+x_pred(4,k)^2);
    L = [m.e_m_(k),m.n_m_(k),m.speed_m_s_(k)]';
    h = [x_pred(1,k),x_pred(2,k),v_]';
    H = [1 0 0 0; 0 1 0 0; 0 0 (x_pred(3,k)/v_) (x_pred(4,k)/v_)];

    Kk = Q_pred(:,:,k) * H'/[R + H * Q_pred(:,:,k) * H'];

    % Step 4: Measurement update
    x_update(:,k) = x_pred(:,k) + Kk * (L - h);

    % Step 5: Covariance update
    Q_update(:,:,k) = (eye(4) - Kk * H) * Q_pred(:,:,k);

    Qx = Q_update(:,:,k);
    x = x_update(:,k);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Smoothing
x_hat = x_update;
Q_hat = Q_update;

for k = N-1:-1:1
    D = Q_update(:,:,k) * Tk' / Q_pred(:,:,k+1);
    x_hat(:,k) = x_update(:,k) + D * (x_hat(:,k+1) - x_pred(:,k+1));
    Q_hat(:,:,k) = Q_update(:,:,k) + D * (Q_hat(:,:,k+1) - Q_pred(:,:,k+1)) * D';
end

D;
x_hat; % diff with true values to know the deviation 
Q_hat;

% After filter and smoothing
error_filt = x_update - tv{:,2:end}';
error_smooth = x_hat - tv{:,2:end}';


% Standard deviations
std_filt = zeros(4, N); % for e, n, ve, vn
std_smooth = zeros(4, N);

for k = 1:N
    std_filt(:,k) = sqrt(diag(Q_update(:,:,k)));
    std_smooth(:,k) = sqrt(diag(Q_hat(:,:,k)));
end

% Plot Measured, Filtered, Smoothed Trajectory
figure;
plot(m.e_m_, m.n_m_, 'k.', 'DisplayName', 'Measured');
hold on;
plot(x_update(1,:), x_update(2,:), 'b-', 'DisplayName', 'Filtered');
plot(x_hat(1,:), x_hat(2,:), 'g-', 'DisplayName', 'Smoothed');
plot(tv{:,2}, tv{:,3}, 'r--', 'DisplayName', 'True');
legend();
xlabel('East [m]'); ylabel('North [m]');
title('Trajectory: Measured vs. Filtered vs. Smoothed vs. True');
grid on;

% Plot Differences (Residuals) to True Values
time = m.Time_s_;

figure;
subplot(2,1,1);
plot(time, error_filt(1,:), 'b-', time, error_smooth(1,:), 'g-', time, zeros(1,N), 'r--');
title('East Position Error'); legend('Filtered', 'Smoothed', 'Zero');
xlabel('Time [s]'); ylabel('Error [m]'); grid on;

subplot(2,1,2);
plot(time, error_filt(2,:), 'b-', time, error_smooth(2,:), 'g-', time, zeros(1,N), 'r--');
title('North Position Error'); legend('Filtered', 'Smoothed', 'Zero');
xlabel('Time [s]'); ylabel('Error [m]'); grid on;

figure;
subplot(2,1,1);
plot(time, error_filt(3,:), 'b-', time, error_smooth(3,:), 'g-', time, zeros(1,N), 'r--');
title('East Velocity Error'); legend('Filtered', 'Smoothed', 'Zero');
xlabel('Time [s]'); ylabel('Error [m/s]'); grid on;

subplot(2,1,2);
plot(time, error_filt(4,:), 'b-', time, error_smooth(4,:), 'g-', time, zeros(1,N), 'r--');
title('North Velocity Error'); legend('Filtered', 'Smoothed', 'Zero');
xlabel('Time [s]'); ylabel('Error [m/s]'); grid on;

% Create table columns
Time = m.Time_s_;
Predicted_E = x_pred(1,:)';
Predicted_N = x_pred(2,:)';
Predicted_Ve = x_pred(3,:)';
Predicted_Vn = x_pred(4,:)';

Filtered_E = x_update(1,:)';
Filtered_N = x_update(2,:)';
Filtered_Ve = x_update(3,:)';
Filtered_Vn = x_update(4,:)';

Smoothed_E = x_hat(1,:)';
Smoothed_N = x_hat(2,:)';
Smoothed_Ve = x_hat(3,:)';
Smoothed_Vn = x_hat(4,:)';

% Create the table
Predicted = table(Time, ...
    Predicted_E, Predicted_N, Predicted_Ve, Predicted_Vn);
Filtered = table(Time, ...
    Filtered_E, Filtered_N, Filtered_Ve, Filtered_Vn);
Smoothed = table(Time, ...
    Smoothed_E, Smoothed_N, Smoothed_Ve, Smoothed_Vn);

% Display the table
disp(Predicted)
disp(Filtered)
disp(Smoothed)





