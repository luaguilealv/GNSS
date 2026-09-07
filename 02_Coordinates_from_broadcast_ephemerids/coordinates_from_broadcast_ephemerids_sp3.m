clc; clear; close all;
format long
%% Data
% Constants
c = 299792458; % speed of light

% Define observation epoch
obs_epoch = datetime(2004, 2, 2, 1, 0, 0);
% Define GPS week start (Sunday 00:00:00 UTC)
gps_week_start = datetime(2004,2,1,0,0,0); 

% Pseudoranges with code P1 and satellites at the specified epoch
P1_pseudoranges = [ 25001256.67943; 23106630.26943; 20664230.28544; 24835459.14243; 23250719.14843; 
    24278147.40443; 21075046.64344; 23040156.09643; 25021425.82543; 24944089.65443; 24444981.72543; 
    20946850.74844];
sat_ids = [24;13;8;21;29;26;10;17;2;28;3;27];
pseudoranges = table(sat_ids,P1_pseudoranges);
pseudoranges = sortrows(pseudoranges,'sat_ids');

% Stellite information according with the rinex document and specification
% of satellites in relatation with the pseudorange
sat_data = readtable("satellite_data.csv", 'Delimiter', ',');
sat_data = sortrows(sat_data,'SatelliteID');
idx = find(ismember(sat_data.SatelliteID,sat_ids));
sat_study = sat_data(idx,:);
sat_number = size(sat_study,1);
% Create the datetime column using existing table columns
sat_study.epoch = datetime(sat_study.Year, sat_study.Month, sat_study.Day, sat_study.Hour, sat_study.Minute, sat_study.Second);

%coordinates provided with precise satellite position 
sat_data_sp3 = readtable("data_igs12561_sp3.csv", 'Delimiter', ',');
idx_sp3 = find(ismember(sat_data_sp3.P,sat_ids));
sat_sp3 = sat_data_sp3(idx_sp3,:);
sat_sp3(:,2:4) = sat_sp3(:,2:4).*1000; % in kilometers

Xs = zeros(sat_number,1);
Ys = zeros(sat_number,1);
Zs = zeros(sat_number,1);

sat_clock_error = zeros(sat_number,1);
%% Compute satellite coordinates and clock correction
for i = 1:sat_number
    % definition of data
    af0 = sat_study.SVClockBias(i);
    af1 = sat_study.SVClockDrift(i);
    af2 = sat_study.SVClockDriftRate(i);
    TGD = sat_study.TGD(i);
    t_oc = mod(seconds(sat_study.epoch(i) - gps_week_start), 604800);
    delta_n = sat_study.Delta_n(i);
    A = (sat_study.sqrtA(i))^2;
    t_oe = sat_study.Toe(i);
    M_0 = sat_study.M0(i);
    e = sat_study.Eccentricity(i);
    omega_0 = sat_study.OMEGA0(i);
    C_us = sat_study.Cus(i);
    C_rs = sat_study.Crs(i);
    C_is = sat_study.Cis(i);
    i_0 = sat_study.i0(i);
    IDOT = sat_study.IDOT(i);
    omega_dot = sat_study.OMEGA_DOT(i);
    C_uc = sat_study.Cuc(i);
    C_rc = sat_study.Crc(i);
    C_ic =sat_study.Cic(i);
    omega = sat_study.omega(i);

    % Step 1: Compute transmission time
    
    % Compute Receiver Time (t_A) in GPS seconds of the week
    t_A_tilda = seconds(obs_epoch - gps_week_start); 

    % Compute Transmission Time (t_s)
    ts_tilda = t_A_tilda - pseudoranges.P1_pseudoranges(i) / c;
     
    % Step 2: Compute satellite clock correction
    dt_sv = af0 + af1 * (ts_tilda - t_oc) + af2 * ((ts_tilda - t_oc)^2);
    dt_s_L1 = dt_sv - TGD;

    % Step 3: Compute system transition time
    ts = t_A_tilda;
    
    % Step 4: Compute eccentric anomaly
    [E_k, M_k,t_k] = step4(A, delta_n, ts, t_oe, M_0, e);

    % Step 5: Compute transmission time
    [d_tr, dt_sv, dt_s_L1] = step5(af0, af1, af2, t_oc, ts_tilda, TGD, e, A, E_k);
    
    
    % Step 6: Repetition of steps 5 and 6 with updated ts
        % Step 4: Compute eccentric anomaly
    [E_k,M_k,t_k] = step4(A, delta_n, ts, t_oe, M_0, e);

        % Step 5: Compute transmission time
    [d_tr, dt_sv, dt_s_L1] = step5(af0, af1, af2, t_oc, ts_tilda, TGD, e, A, E_k);
    
    % Step 7: Compute satellite coordinates
    [xk, yk, zk] = step7(A, e, E_k,omega, omega_0, C_uc, C_rc, C_ic, C_us, C_rs, C_is, i_0, IDOT, t_k, omega_dot,t_oe);

    Xs(i) = xk;
    Ys(i) = yk;
    Zs(i) = zk;

    sat_clock_error(i) = dt_s_L1;
end

res = table (sat_study.SatelliteID, Xs, Ys, Zs, sat_clock_error);
disp(res);
%% Difference with precise satellite positions
% Compute error
error_X = res.Var2 - sat_sp3.X;
error_Y = res.Var3 - sat_sp3.Y;
error_Z = res.Var4 - sat_sp3.Z;
errors = table(error_X,error_Y,error_Z);
disp(errors);
%% Definition of the functions
% Function to compute Mean anomaly and Kepler's Equation for eccentric anomaly
function [E_k, M_k,t_k] = step4(A, delta_n, ts, t_oe, M_0, e)
    mu = 3.986005e14; % Earth's universal gravitational parameter 

    % Mean motion
    n_0 = sqrt(mu / A^3);
    % Correct mean motion
    n = n_0 + delta_n;

    % Time from ephemeris reference epoch
    t_k = definetk(ts,t_oe);
        
    % Mean anomaly
    M_k = M_0 + n * t_k;
    
    % Iterative solution Kepler's Equation for eccentric anomaly
    E_k = M_k;
    for iter = 1:5 
        E_0 = E_k;
        E_k = M_k + e * sin(E_0);
        if abs(E_k - E_0) < 1e-13
            break;
        end
    end
end

% Function to define Time from ephemeris reference epoch
function tk = definetk(ts, t_oe)
    tk = (ts - t_oe);
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -302400
        tk = tk + 604800;
    end
end

% Function to compute satellite clock correction
function [d_tr, dt_sv, dt_s_LI] = step5(af0, af1, af2, t_oc, ts_tilda, T_GD, e, A, E_k)
    F = -4.442807633e-10; % Relativistic correction constant
    d_tr = F * e * sqrt(A) * sin(E_k);
    dt_sv = af0 + af1 * (ts_tilda - t_oc) + af2 * (ts_tilda - t_oc)^2 + d_tr;
    dt_s_LI = dt_sv - T_GD;
end

% Function to compute satellite coordinates
function [xk, yk, zk] = step7(A, e, E_k,omega, omega_0, C_uc, C_rc, C_ic, C_us, C_rs, C_is, i_0, IDOT, t_k, omega_dot,t_oe)
    Omega_dot_e = 7.2921151467e-5; % Earth's rotation rate

    % True anomaly
    sin_v_k = sqrt(1 - e^2) * sin(E_k) / (1 - e * cos(E_k));
    cos_v_k = (cos(E_k) - e) / (1 - e * cos(E_k));
    v_k = atan2(sin_v_k, cos_v_k);

    % Argument of latitude
    phi_k = v_k + omega;

    % Latitude, Radius, Inclination correction 
    delta_u_k = C_us * sin(2 * phi_k) + C_uc * cos(2 * phi_k);
    delta_r_k = C_rs * sin(2 * phi_k) + C_rc * cos(2 * phi_k);
    delta_i_k = C_is * sin(2 * phi_k) + C_ic * cos(2 * phi_k);

    % Corrected argument of latitude, radius, inclination
    u_k = phi_k + delta_u_k;
    r_k = A * (1 - e * cos(E_k)) + delta_r_k;
    i_k = i_0 + delta_i_k + IDOT * t_k;

    % Positions in orbital plane
    x_prime_k = r_k * cos(u_k);
    y_prime_k = r_k * sin(u_k);

     % Corrected longitude of ascending node
    omega_k = omega_0 + (omega_dot - Omega_dot_e) * t_k - Omega_dot_e * t_oe;

    % Cartesian coordinates
    xk = x_prime_k * cos(omega_k) - y_prime_k * cos(i_k) * sin(omega_k);
    yk = x_prime_k * sin(omega_k) + y_prime_k * cos(i_k) * cos(omega_k);
    zk = y_prime_k * sin(i_k);
end 