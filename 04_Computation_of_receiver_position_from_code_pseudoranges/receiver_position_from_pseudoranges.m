clc; clear; close all;
format long
% read satellites coordinates
satellites = xlsread("testdata.xlsx");
sat_number = size(satellites,1);
dt_s_L1 = satellites(:,5);

% Define observation epoch
obs_epoch = datetime(2004, 2, 2, 1, 0, 0);
% Define GPS week start (Sunday 00:00:00 UTC)
gps_week_start = datetime(2004,2,1,0,0,0);

% Pseudoranges with code P1 and satellites at the specified epoch
P1_pseudoranges = [ 25021425.82543; 24444981.72543; 20664230.28544; 21075046.64344; 23106630.26943; 
    23040156.09643; 24835459.14243; 25001256.67943; 24278147.40443; 20946850.74844; 24944089.65443; 
    23250719.14843];
    
% Step 8: compute tropospheric correction
% Tropospheric delay in meters for each satellite
Tropospheric_delay = [9.25963303386247; 10.8638032902129; 2.65787939593671; 2.81147015429228; 5.02959675957772;
    5.59603189513188; 13.2116073589126; 14.7475814061828; 12.9837949367863; 2.66518532200263; 13.0073814676033;
    5.1778396245971];

% Step 9: compute ionospheric correction
% Ionospheric delay in meters for each satellite
Ionospheric_delay = [1.012141327186487; 1.069019624266762; 0.442846108799202; 0.465181639969468;
   0.733232317847840; 0.785925529559771; 1.126368176698085; 1.154557609247377; 1.121732297962704; 
   0.443919134357337; 1.122220642853633; 0.747589301888461];

% Approximate coordinates [X; Y; Z]
x_A0 = [3104219; 998383; 5463290];


% Initial values
Omega_dot_e = 7.2921151467e-5; % Earth's rotation rate
theta = Omega_dot_e*(satellites(:,6) - seconds(obs_epoch - gps_week_start));
Xs_prime = satellites(:,2).*cos(theta)-satellites(:,3).*sin(theta);
Ys_prime = satellites(:,2).*sin(theta)+satellites(:,3).*cos(theta);
Zs_prime = satellites(:,4);

% Algorithm to compute receiver's position for the other steps
epsilon = 1e-5;
max_iter = 200; % safety to avoid infinite loops
iter = 0;
converged = false;
v_0 =  zeros(size(P1_pseudoranges));  % initial residuals;  % initial residuals
x_A = x_A0;  % start with initial approx. coordinates

while ~converged && iter < max_iter
    iter = iter + 1;
    rho_a0_s = step10(Xs_prime, Ys_prime, Zs_prime, x_A(1), x_A(2), x_A(3));% Step 10
    L = step11(rho_a0_s, dt_s_L1, Ionospheric_delay, Tropospheric_delay, P1_pseudoranges);% Step 11  
    [A, ax_s, ay_s, az_s] = step12(Xs_prime, Ys_prime, Zs_prime, x_A(1), x_A(2), x_A(3), rho_a0_s);% Step 12
    [Q, X_hat, v] = step13(A, L);% Step 13
    x_A = step14(x_A, X_hat);% Step 14
    % Step 15: check convergence
    if abs(v'*v - v_0'*v_0) < epsilon
        converged = true;
    end
    % Update v_0 for next iteration
    v_0 = v;
end

s0 = sqrt(v'*v/(sat_number-4));
sX = s0*sqrt(Q(1,1));
sY = s0*sqrt(Q(2,2));
sZ = s0*sqrt(Q(3,3));
sc = s0*sqrt(Q(4,4));
PDOP = sqrt(Q(1,1)+Q(2,2)+Q(3,3));

fprintf("Estimated receiver position:\nX = %.3f m\nY = %.3f m\nZ = %.3f m\n", x_A(1), x_A(2), x_A(3));
fprintf("Receiver clock bias = %.3f m \n", X_hat(4));
fprintf("PDOP = %.3f\n", PDOP);

function dist = step10(Xs,Ys,Zs,Xa0,Ya0,Za0)
% Step 10: compute approximate distance
    dist = sqrt((Xs-Xa0).^2+(Ys-Ya0).^2+(Zs-Za0).^2);
end

function L = step11(rho_a0_s,dt_s_L1,I_a,T_a,P_a)
% Step 11: compute vector L
    c = 299792458; % speed of light
    L = P_a - rho_a0_s + c*dt_s_L1 - I_a - T_a;
end

function [A,ax_s,ay_s,az_s] = step12(Xs,Ys,Zs,Xa0,Ya0,Za0,rho_a0_s)
% Step 12: compute elements of matrix A
    ax_s = -(Xs-Xa0)./rho_a0_s;
    ay_s = -(Ys-Ya0)./rho_a0_s;
    az_s = -(Zs-Za0)./rho_a0_s;
    uno = ones(size(ax_s));
    A = [ax_s ay_s az_s uno];
end

function [Q,X_hat,v] = step13(A,L)
% Step 13: compute unknown parameters
    Q = (A'*A) \ eye(size(A,2));
    X_hat = Q*A'*L;
    v = L - A*X_hat;
end

function [x_A] = step14(x0,X_hat)
% Step 14: update receivers coordinates
    x_A = x0 - X_hat(1:3,:);
end