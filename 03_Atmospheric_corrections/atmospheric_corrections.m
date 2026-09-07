clc; clear; close all;
format long
% read satellites coordinates
satellites = xlsread("testdata.xlsx");
sat_number = size(satellites,1);

% define the station coordinates in cartesian and elliptic
X_station = 3104219.4530;
Y_station = 998383.982;
Z_station = 5463290.5080;
% GRS80 ellipsoid
a = 6378137.0;
f = 1 / 298.257222101;
[lat_station, lon_station ,~] = cartesiantogeodetic(X_station, Y_station, Z_station, a, f);
h_station = 62 + 23;

% compute satellites azimuths
[~, ~, ~, zeniths, ~] = enucoordinates(satellites, lat_station, lon_station, X_station,Y_station,Z_station);

% Constants
T0_C = 18;  % Temperature in Celsius
P = 1013;   % Pressure in mbar
RH = 0.5;   % Relative Humidity (50%)

% Convert temperature to Kelvin
T0_K = T0_C + 273.15;

% Compute Partial Pressure of Water Vapor (e) in hPa = mbar
es = RH * 6.108 * exp((17.15 * T0_K - 4684) / (T0_K - 38.45));

%% Saastamoinen model
% Without correction terms
delta_trop_uc = arrayfun(@(z) troposphericdelayuncorrected(P , T0_K, es, z), zeniths);

delta_trop_uc_test = troposphericdelayuncorrected(1013.25 , 273.16, 0, 0); % check tropospheric delay

% With correction terms
delta_trop = arrayfun(@(z) troposphericdelay(P ,T0_K, es, z, h_station), zeniths);

% check tropospheric delay
delta_trop_test = troposphericdelay(1013.25 , 273.16, 0, 0, 0);

tab = table(satellites(:,1), delta_trop_uc, delta_trop);
disp(tab);


%% Change in temperature and humidity

temp = linspace(-30, 30, 100) + 273.15;
hum = linspace(0, 1, 100);
ztd_temp = zeros(length(temp), sat_number);
ztd_hum = zeros(length(hum), sat_number);


% only temperature
for j = 1:sat_number
    for i = 1:length(temp)
         es_t = RH * 6.108 * exp((17.15 * temp(i) - 4684) / (temp(i) - 38.45));
         ztd_temp(i,j) = troposphericdelay(P , temp(i), es_t, zeniths(j), h_station);
    end
end

% Plot results
figure;
for j = 1:sat_number
    subplot(2, ceil(sat_number / 2), j);
    plot(temp, ztd_temp(:, j), 'b', 'LineWidth', 1.5);
    xlabel('Temperature (K)');
    ylabel('ZTD (m)');
    title(sprintf('Zenith Angle %.1f°', zeniths(j)));
    grid on;
end
sgtitle('Effect of Temperature on Zenith Tropospheric Delay'); % Global title

% only humidity
for j = 1:sat_number
    for i = 1:length(hum)
         es_h = hum(i) * 6.108 * exp((17.15 * T0_K - 4684) / (T0_K - 38.45));
         ztd_hum(i,j) = troposphericdelay(P , T0_K, es_h, zeniths(j), h_station);
    end
end
% Plot results
figure;
for j = 1:sat_number
    subplot(2, ceil(sat_number / 2), j);
    plot(hum, ztd_hum(:, j), 'b', 'LineWidth', 1.5);
    xlabel('Humidity (%)');
    ylabel('ZTD (m)');
    title(sprintf('Satellite %.1f°', zeniths(j)));
    grid on;
end
sgtitle('Effect of Humidity on Zenith Tropospheric Delay'); % Global title


% Change in temperature and humidity
ZTD = zeros(length(temp), length(hum));

% Compute ZTD for each Temperature & Humidity Combination
for k = 1:sat_number
    for i = 1:length(temp)
        for j = 1:length(hum)
            RH = hum(j);  % Current Relative Humidity
            % Compute water vapor partial pressure (es_t)
            es_c = RH * 6.108 * exp((17.15 * temp(i) - 4684) / (temp(i) - 38.45));
            % Compute Zenith Tropospheric Delay (ZTD)
            ZTD(i, j, k) = troposphericdelay(P, temp(i), es_c, zeniths(k), h_station);
        end
    end
end

% Plot Heatmaps for Each Zenith Angle
figure;
for k = 1:sat_number
    subplot(3, 4, k); % Create 3x4 grid for 12 plots
    imagesc(hum * 100, temp, ZTD(:,:,k)); % Convert humidity to percentage
    colorbar;
    xlabel('Humidity (%)');
    ylabel('Temperature (K)');
    title(sprintf('Zenith Angle %.1f°', zeniths(k)));
    set(gca, 'YDir', 'normal'); % Ensure temperature is plotted correctly
end
sgtitle('Effect of Temperature & Humidity on Zenith Tropospheric Delay'); % Global title

%% Ionospheric delay
f_L1 = 1575420000; %Hz, frequency  
h_I = 350000; % m, mean ionospheric height
R_E = 6371000; % m, Earth mean radius
TEC = 2.5e16; % el/m2, Total electron content

I_delay = arrayfun(@(z) ionosphericdelay(f_L1, TEC, h_I, R_E, z), zeniths);

disp(I_delay);

%% Define the functions

% Function to compute Saastamoinen Model without correction
function delta_trop = troposphericdelayuncorrected(p , T, e, z)
    delta_trop = (0.002277/cosd(z))*(p + e * (0.05 + (1255 / T)) - (tand(z))^2);
end
% Function to compute Saastamoinen Model with correction
function delta_trop = troposphericdelay(p, T, e, z, station_height)
    height_levels = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]; % in kilometers
    height_levels_R = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]; % in kilometers
    % correction term B from Table 5.4
    B_values = [1.156, 1.079, 1.006, 0.938, 0.874, 0.813, 0.757, 0.654, 0.563];
    % correction term δR from Table 5.5 
    R_matrix = [
        0.003, 0.003, 0.002, 0.002, 0.002, 0.002, 0.001, 0.001;
        0.006, 0.006, 0.005, 0.005, 0.004, 0.003, 0.003, 0.002;
        0.012, 0.011, 0.010, 0.009, 0.008, 0.006, 0.005, 0.004;
        0.020, 0.018, 0.017, 0.015, 0.013, 0.011, 0.009, 0.007;
        0.031, 0.028, 0.025, 0.023, 0.021, 0.017, 0.014, 0.011;
        0.039, 0.035, 0.032, 0.029, 0.026, 0.021, 0.017, 0.014;
        0.050, 0.045, 0.041, 0.037, 0.033, 0.027, 0.022, 0.018;
        0.065, 0.059, 0.054, 0.049, 0.044, 0.036, 0.030, 0.024;
        0.075, 0.068, 0.062, 0.056, 0.051, 0.042, 0.034, 0.028;
        0.087, 0.079, 0.072, 0.065, 0.059, 0.049, 0.040, 0.033;
        0.102, 0.093, 0.085, 0.077, 0.070, 0.058, 0.047, 0.039;
        0.111, 0.101, 0.092, 0.083, 0.076, 0.063, 0.052, 0.043;
        0.121, 0.110, 0.100, 0.091, 0.083, 0.068, 0.056, 0.047;
    ];

    zenith_angles = [60, 66, 70, 73, 75, 76, 77, 78, 78.3, 79, 79.3, 79.45, 80];

    % Create interpolation grid
    [H, Z] = meshgrid(height_levels_R, zenith_angles);

    % Interpolation
    B = interp1(height_levels, B_values, station_height / 1000, 'linear', 'extrap'); 
    delta_R = interp2(H, Z, R_matrix, station_height / 1000, z, 'linear');

    % Ensure delta_R is within reasonable bounds
    if z < 60
        delta_R = 0.001;
    elseif z < 30
        delta_R = 0;
    elseif z > 80
       if station_height < 1.5
          delta_R = 0.1055;
       else
           delta_R = 0.0635;
       end
    end   

    % Refined model
    delta_trop = (0.002277 / cosd(z)) * (p + (0.05 + 1255 / T) * e - B * (tand(z))^2) + delta_R;
end


% function to compute Ionospheric delay
function I_f = ionosphericdelay(f, TEC, h, R, z)
    OF = (1 - ((R * sind(z))/(R + h))^2)^(-1/2);
    I_f_z = TEC * 40.3 / (f^2);
    I_f = I_f_z * OF;
end


% define enu coordinates, azimuths, ...
function [enu, azimuths, slant_distances, zenith_angles, elevation_angles] = enucoordinates(satellites, lat_given, lon_given, X,Y,Z)
    % funtion to compute azimuts, distances, zenith and elevation angles
    % transforme coordinates to radians
    lat = deg2rad(lat_given);
    lon = deg2rad(lon_given);

    % number of observations
    num_sats = size(satellites, 1);

    % from ECEF to ENU coordinates
    % Rotation matrix 
    % R = Rx(90º-lat)Rz(lon+90º)
    R = [-sin(lon), cos(lon), 0;
         -sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat);
          cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)];

    % ENU coordinates for each satellite
    % ENU = R*CTRS; CTRS = SV - a; [X_SV - Xa; Y_SV - Ya; Z_SV - Za]

    enu = zeros(num_sats, 3);
    azimuths = zeros(num_sats, 1);
    slant_distances = zeros(num_sats, 1);
    zenith_angles = zeros(num_sats, 1);
    elevation_angles = zeros(num_sats, 1);

    for i = 1 : num_sats
        % Satellite 
        X_SV = satellites(i, 2);
        Y_SV = satellites(i, 3);
        Z_SV = satellites(i, 4);

        enu_SV = R * [X_SV - X; Y_SV - Y; Z_SV - Z];
        enu(i, :) = enu_SV';
        e = enu_SV(1);
        n = enu_SV(2);
        u = enu_SV(3);

        % azimuth, slant distance, zenith angle, elevation angle
        azimuths(i) = mod(rad2deg(atan2(e, n)), 360);
        slant_distances(i) = sqrt(n^2 + e^2 + u^2);
        zenith_angles(i) = acosd(u/slant_distances(i));
        elevation_angles(i) = 90 - zenith_angles(i);
    end
end

% function to convert cartesian coordinates to geodetic
function [lat, lon , h] = cartesiantogeodetic(X, Y, Z, a, f)
% function to convert cartesian coordinates to geodetic
    % f = (a - b)/a
    % semi minor axis
    b = a - f*a;

    % first numerical eccentricity
    e2 = (a^2 - b^2)/(a^2);

    % longitude
    lon = atan2d(Y, X);

    % iterative solution
    % radious of a parallel
    p = sqrt(X^2 + Y^2);
    % latitude
    lat0 = atan2d(Z, p * (1 - e2));

    % latitude and height
    lat = 0;
    while abs(lat - lat0) > 1e-12  % Convergence criterion
        N = a^2 / sqrt( a^2*cos(lat0)^2 + b^2*sin(lat0)^2);
        h = p / cosd(lat0) - N;
        lat0 = atan2d(Z, p * (1 - e2 * N / (N + h)));
        lat = lat0;
    end   
end