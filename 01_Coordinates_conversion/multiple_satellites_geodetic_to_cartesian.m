clear; clc; close all;
%% Data
% GRS80 ellipsoid
a = 6378137.0;
f = 1 / 298.257222101;

% geodetic coordinates
lat_given = dms2degrees([59 20 59]);
lon_given = dms2degrees([18 4 10]);

% ellipsoidal height h
% geoidal height N
% orthometric height, height above sea level
H = 62;
N = 23;
h = H + N; % h

% Convert geodetic to Cartesian coordinates
[X,Y,Z] = geodetictocartesian(lat_given, lon_given, h, a, f);

% Load Satellite Positions in kilometers 
% columns 2, 3, 4 of the table are coordinates, first column identification of satellite 
satellites = readmatrix("GPS satellite positions from sp3 file.xlsx");
% Convert km to meters
satellites(:, 2:4) = satellites(:, 2:4) * 1e3;

%% control point
% test position (100 meters above the origin position)
h_test = h + 100; % New height
[Xt, Yt, Zt] = geodetictocartesian(lat_given, lon_given, h_test, a, f);

% ENU coordinates of test point relative to the origin
enu2 = enucoordinates([1,Xt, Yt, Zt], lat_given, lon_given, X,Y,Z);

%% results
[enu, azimuths, slant_distances, zenith_angles, elevation_angles] = enucoordinates(satellites, lat_given, lon_given, X,Y,Z);

%sumarize data
satellite_number = satellites(:,1);
sat_info = table(satellite_number,enu,elevation_angles,azimuths,slant_distances);

% visible satellites
v_idx = find(elevation_angles > 0);
nv_idx = find(elevation_angles <= 0);
visible_sats = length(v_idx);

disp('Number of visible satellites:');
disp(visible_sats);
disp('Visible satellites:');
disp(sat_info(v_idx,:));

%% Visualizations
% satellites position in 2D (skyplot(azimut,elevation))
figure(1);
skyplot(azimuths(v_idx), elevation_angles(v_idx)); 
title('Visualization of the visible satellites');


% 3D Plot
figure(2);
[Xs, Ys, Zs] = sphere(50); 
earth_surface = surf(Xs * 6371e3, Ys * 6371e3, Zs * 6371e3);
set(earth_surface, 'EdgeColor', 'none');
hold on 
% plot3(enu(v_idx, 1), enu(v_idx, 2), enu(v_idx, 3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% plot3(enu(nv_idx, 1), enu(nv_idx, 2), enu(nv_idx, 3), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
plot3(satellites(v_idx, 2), satellites(v_idx, 3), satellites(v_idx, 4), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot3(satellites(nv_idx, 2), satellites(nv_idx, 3), satellites(nv_idx, 4), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
plot3(X,Y,Z, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Visualization of satellites');
legend('earth','visible satellites','nonvisible satellites','observer')
grid on; axis equal;
hold off;

%% functions
function [X,Y,Z] = geodetictocartesian(lat, lon, h, a, f)
% function to convert geodetic coordinates to cartesian coordinates
    % lat = phy, lon = lambda
    % Convert latitude and longitude from degrees to radians
    lat = deg2rad(lat);
    lon = deg2rad(lon);

    % f = (a - b)/a
    % semi minor axis
    b = a - f*a;

    % first numerical eccentricity
    e2 = (a^2 - b^2)/(a^2);

    % radius of curvature in the prime vertical
    N = a^2 / sqrt( a^2*cos(lat)^2 + b^2*sin(lat)^2);

    % Cartesian coordinates
    X = (N + h) * cos(lat) * cos(lon);
    Y = (N + h) * cos(lat) * sin(lon);
    Z = (N + h - e2*N) * sin(lat);
end

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