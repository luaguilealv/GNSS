clear; clc; close all;
% GRS80 ellipsoid
a = 6378137.0;
f = 1 / 298.257222101;

% geodetic coordinates
lat_ref = dms2degrees([59 20 59]);
lon_ref = dms2degrees([18 4 10]);

% ellipsoidal height h
% geoidal height N
% orthometric height, height above sea level
H = 62;
N = 23;
h_ref = H + N; % h

% Cartesian coordinates
X_ref = 3098917.24;
Y_ref = 1011053.41;
Z_ref = 5463972.35;

% Convert geodetic to Cartesian
[X,Y,Z] = geodetictocartesian(lat_ref, lon_ref, h_ref, a, f);
disp('Computed Cartesian coordinates:');
disp([X, Y, Z]);

% difference between given coordinates and computed
diff_X = X - X_ref;
diff_Y = Y - Y_ref;
diff_Z = Z - Z_ref;

disp('Differences in Cartesian coordinates:');
disp([diff_X, diff_Y, diff_Z]);

% Convert Cartesian to geodetic
[lat, lon , h] = cartesiantogeodetic(X_ref, Y_ref, Z_ref, a, f);
disp('Computed geodetic coordinates:');
disp([lat, lon, h]);

% Reconvert Cartesian to geodetic from the computed coordinates
[lat_com, lon_com , h_com] = cartesiantogeodetic(X, Y, Z, a, f);
disp('Reconverted geodetic coordinates:');
disp([lat_com, lon_com, h_com]);

% difference between given coordinates and computed
diff_lat = lat_ref - lat_com;
diff_lon = lon_ref - lon_com;
diff_h = h_ref - h_com;

disp('Differences in geodetic coordinates:');
disp([diff_lat, diff_lon, diff_h]);

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