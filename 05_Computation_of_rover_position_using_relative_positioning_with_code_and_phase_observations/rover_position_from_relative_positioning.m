clc; clear; close all;
format long

% read data
REF = readtable('Reference.xlsx');
ROV = readtable('Rover.xlsx');

lambda = 0.190293673; % m

% Receivers coordinates [X,Y,Z] [m]
ref_approx = [3098896.003, 1011042.18, 5463974.974];
rov_approx = [3101687.4, 1012242.4, 5462199.3];

% Step 1: Synchronize observables from both receivers using equation (15)
P_REF = REF.P1_m_ + REF.Rho_dot.*REF.ReceiverClockError_m_/299792458;
phi_REF = lambda*REF.L1_cycles_ + REF.Rho_dot.*REF.ReceiverClockError_m_/299792458;

P_ROV = ROV.P1_m_ + ROV.Rho_dot.*ROV.ReceiverClockError_m_/299792458;
phi_ROV = lambda*ROV.L1_cycles_ + ROV.Rho_dot.*ROV.ReceiverClockError_m_/299792458;

% Step 2: Compute single differences for both code and phase using equation (16)
phi_SD = phi_ROV - phi_REF;
P_SD = P_ROV - P_REF;

% Step 3: Compute double differences using (19) and use satellite 1 as the reference satellite
[n,m] = size(P_SD);
e = ones(n-1,m);
phi_DD = phi_SD(1)*e - phi_SD(2:end);
P_DD = P_SD(1)*e - P_SD(2:end);

% % Step 4: Compute the coefficients aX, aY, aZ using equation (13) and the
% % single difference 𝜌𝐴𝐵,0𝑠 by subtracting the two sets of distances between receivers and satellites.
v = ones(size(P_SD));
aX = -(ROV.X_m_ - rov_approx(1)*v)./ROV.TopocentricDistanceRho;
aY = -(ROV.Y_m_ - rov_approx(2)*v)./ROV.TopocentricDistanceRho;
aZ = -(ROV.Z_m_ - rov_approx(3)*v)./ROV.TopocentricDistanceRho;

rho_SD = ROV.TopocentricDistanceRho-REF.TopocentricDistanceRho;

% Step 5: Compute the double differenced coefficients using equation (23) and 𝜌𝐴𝐵,0 𝑠𝑡 using equation (20).
aX_DD = aX(1)*e - aX(2:end);
aY_DD = aY(1)*e - aY(2:end);
aZ_DD = aZ(1)*e - aZ(2:end);


rho_DD = rho_SD(1)*e - rho_SD(2:end); 

% Step 6: Fill in matrix A and the vector L and compute a least squares solution of equation (22) using equation (27). 
% Compute A
A_block = [aX_DD, aY_DD, aZ_DD];
I = lambda*eye(n-1);
Z = zeros(size(I));
A = [A_block I; A_block Z];

% Compute L
L = [phi_DD-rho_DD;P_DD-rho_DD];

% Compute P
sigma_phi = 0.002; % phase std in m
sigma_P = 0.3;  % code std in m
n_dd = length(e);
M = -1*ones(n_dd) + eye(n_dd) + n_dd*eye(n_dd);
P_P = (1/(2*sigma_phi^2*(n_dd+1)))*M;
P_phi = (1/(2*sigma_P^2*(n_dd+1)))*M;
P = blkdiag(P_P, P_phi);

% Compute X
X = (A' * P * A) \ (A' * P * L);

% Step 7: Correct the approximate position of the rover with the corrections for X, Y, and Z given in the x-vector.
% New coordinates
Xrov_new = rov_approx(1) + X(1);
Yrov_new = rov_approx(2) + X(2);
Zrov_new = rov_approx(3) + X(3);

% Step 8: Compute the variance-covariance matrix (Qahat) from the least squares estimation of the float 
% ambiguities (ahat). Then extract the ndd x ndd sub-matrix for the ambiguities. 
% ndd = number of double differenced ambiguities
Qx = inv(A' * P * A);  % Full covariance matrix (size: n_x × n_x)

idx_N = 4:(3+n_dd);  

Qahat = Qx(idx_N, idx_N);  % n_dd x ndd matrix
ahat = X(idx_N); % float ambiguity estimates

% Step 9: Call the LAMBDA function 
[afixed, sqnorm, Ps, Qzhat, Z, nfixed, mu] = LAMBDA(ahat, Qahat);

% Step 10: Fill in the L vector using only the phase double differences and subtract the fixed ambiguities
L_fixed = (phi_DD - rho_DD) - afixed(:,1)*lambda; % best set of ambiguities
L_next = (phi_DD - rho_DD) - afixed(:,2)*lambda;  % next best set of ambiguities

% Step 11: Fill in the A matrix only for the phase double differences
A_fixed = A(1:n_dd,1:3);

% Step 12: Fill in the P matrix only for the phase double differences, using the same weights as before
P_fixed = P_phi;

% Step 13: Solve for X and correct the approximate position for ROV as it was given to estimate the
% position for ROV with the fixed ambiguities.
x_fixed = (A_fixed' * P_fixed * A_fixed) \ (A_fixed' * P_fixed * L_fixed); % best set of ambiguities
x_next = (A_fixed' * P_fixed * A_fixed) \ (A_fixed' * P_fixed * L_next); % next best set of ambiguities

XR_final = rov_approx(1) + x_fixed(1);
YR_final = rov_approx(2) + x_fixed(2);
ZR_final = rov_approx(3) + x_fixed(3);

% Step 14: Perform ambiguity validation. 
% Above call to the LAMBDA function the best and the next best set of ambiguities are returned (afixed).
% Compute the residuals, (v), with both the best and the next best solution of fixed ambiguities. 
% Then calculate the ratio test variable and evaluate the set of fixed ambiguities.
v1 = L_fixed - A_fixed * x_fixed;
v2 = L_next - A_fixed * x_next;

ratio = (v2'* P_fixed * v2) / (v1'* P_fixed * v1);

if ratio > 3  % Typically
    disp('Accept');
    disp(ratio);
else
    disp('Reject');
    disp(ratio);
end

% Display float vs fixed ambiguities
fprintf('\nFloat vs Fixed Ambiguities:\n');
fprintf('  #   Float Ambiguity    Fixed Ambiguity    Difference\n');
for i = 1:n_dd
    fprintf('%2d   %14.6f   %14d   %10.6f\n', i, ahat(i), afixed(i,1), ahat(i) - afixed(i,1));
end

% Extract standard deviations of the final coordinates (from Qx)
Qx_fixed = inv(A_fixed' * P_fixed * A_fixed);  % 3x3 covariance matrix
std_devs = sqrt(diag(Qzhat));  % [σ_X; σ_Y; σ_Z]

% Display final coordinates and std devs
fprintf('\nFinal Rover Coordinates (Fixed Ambiguities):\n');
fprintf('X: %.4f ± %.4f m\n', XR_final, std_devs(1));
fprintf('Y: %.4f ± %.4f m\n', YR_final, std_devs(2));
fprintf('Z: %.4f ± %.4f m\n', ZR_final, std_devs(3));
