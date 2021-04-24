%% Read data
% The input-output data was collected for scaled model of an aircraft kept
% in a wind tunnel. The file roll_data_wind_tunnel_test.csv contains the
% control input commands and output measurements.
Data = readmatrix('roll_data_wind_tunnel_test.csv');

% Column 1: Roll-angle measurements (radians).
phi = Data(:,1);

% Column 2: Roll-rate measurements (radians/second).
p = Data(:,2);

% Column 3: Aileron deflection angle (radians).
deltaA = Data(:,3);

% We captured 1500 samples at 100Hz sample rate
Nsamps = 1500;
fsamp = 100;
Tsamp = 1/fsamp;
t = (0:Nsamps-1)*Tsamp;

%% System-Identification
% Assume the first order model from aeilron deflection to the roll-rate to
% estimate model parameters. See attached note for more details. The model
% to identfy has a transfer function G(s) = k/(s+a). Select the approach to
% perform system identification
approach_case = 1;
switch approach_case
    case 1
        % Use system identification app to obtain the following
        tfdata  = iddata(p,deltaA,Tsamp);
        Gest = tfest(tfdata,1,0);
        a = Gest.Denominator(2);
        k = Gest.Numerator;
        
    case 2
        % Use Euler approximation to write discrete-time dynamics and
        % estimate coefficients using least-square
        Z = [p(1:end-1), deltaA(1:end-1)];
        b = p(2:end);
        pStar = Z\b;
        a = (1-pStar(1))/Tsamp;
        k = pStar(2)/Tsamp;
end

% State-space matrices: states are roll-angle (phi) and roll-rates (p)
A = [0 1; 0 -a];
B = [0; k];
C = eye(2);
D = 0;
G = ss(A,B,C,D);

%% Linear Simulation with Identified Model
% lsim
[xs,ts] = lsim(G,deltaA,t);

% Measured data
xm = [phi,p];

%% Plot Data
figure(1),
plot(t,xm(:,1),'b',ts,xs(:,1),'c','LineWidth',2);grid on;box on;
xlabel('Time (sec)','FontSize',12);ylabel('Roll Angle \phi (rad)','FontSize',12);
legend('Measured','Computed','FontSize',12);

figure(2),plot(t,xm(:,2),'r',ts,xs(:,2),'k','LineWidth',2);grid on;box on;
xlabel('Time (sec)','FontSize',12);ylabel('Roll Rate p (rad/sec)','FontSize',12);
legend('Measured','Computed','FontSize',12);

figure(3),plot(t,deltaA,'g','LineWidth',2);grid on;box on;
xlabel('Time (sec)','FontSize',12);ylabel('\delta_A (rad)','FontSize',12);

%% Compute MSE
% Compute suqared mean error in roll angle estimate
% Note that e1RMS is less than specified number 5e-4
e1RMS = sum((xm(:,1) - xs(:,1)).^2)/Nsamps;

% Compute suqared mean error in roll-rate estimate
e2RMS = sum((xm(:,2) - xs(:,2)).^2)/Nsamps;