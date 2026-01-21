function [PVx, PVy] = kalmanFilter(curDetX, curDetY, curDetZ, dt_)
%KALMANFILTER Calculates Primary vertex using kalman filter

% State is [x; vx; y; vy; z; vz]
totalPions = size(curDetX,1);
PVx = zeros(totalPions,1);
PVy = zeros(totalPions,1);
for i = 1:totalPions
    xPos = curDetX(i,:);
    yPos = curDetY(i,:);
    zPos = curDetZ;
    time = dt_(i,:);

    % State transition matrix A = Block diagonal matrix
    % with the [1 dt; 0 1] block repeated for the 
    % x, y, and z spatial dimensions.
    dt = (time(2)+time(1))/2;
    A1D = [1 dt; 0 1];
    A = kron(eye(3),A1D); % State transiton model
    
    % Initial estimate uncertainty (Covariance matrix)
    % Assume intial guess is not great 1 mm error
    P0 = eye(6)*0.0000001;
    
    % Process noise matrix for a discrete model
    q=0.001;
    Q1D = [(dt^4)/4, (dt^3)/2;
        (dt^3)/2, dt^2];
    Q = kron(eye(3),Q1D)*q^2;
    
    % Observation matrix
    % z_n = H*x_n 
    % z measures x,y,z position so dimension 3x1
    % x_n has dimension 6x1 so H is 3x6
    H = [1, 0, 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0;
        0, 0, 0, 0, 1, 0];
    
    % Measurement uncertainty (covariance) matrix
    % SciFi have resolution of 100 um 
    sigma = 0.01;
    R = sigma^2*eye(3);

    % Don't know particle info, so put in reasonable guess
    x0 = [0;0;0;0;0;2.55e8];

    % Predict next states 
    x1 = A*x0;
    P1 = A*P0*A' + Q;

    % First iteration
    z1 = [xPos(1);yPos(1);zPos(1)];     % Measured values
    K1 = P1*H'*(H*P1*H' + R)^-1;        % Kalman Gain calculation
    x1_ = x1 + K1*(z1-H*x1);            % Estimate current state
    P1_ = (eye(6)-K1*H)*P1*(eye(6) -K1*H)' + K1*R*K1';  % Update the estimate uncertainty
    x2 = A*x1_;                         % Predict next state
    P2 = A*P1_*A' + Q;                  % Predict next state
    
    % Iteration 2
    z2 = [xPos(2);yPos(2);zPos(2)];     % Next measured values
    K2 = P2*H'*(H*P2*H' + R)^-1;        % Kalman gain calculation
    x2_ = x2 + K2*(z2-H*x2);            % Estimate current state
    P2_ = (eye(6)-K2*H)*P2*(eye(6) -K2*H)' + K2*R*K2';  % Update the estimate uncertainty
    x3 = A*x2_;                         % Predict next state
    P3 = A*P2_*A' + Q;                  % Predict next state

    % Iteration 3
    z3 = [xPos(3);yPos(3);zPos(3)];     % Next measured values
    K3 = P3*H'*(H*P3*H' + R)^-1;        % Kalman gain calculation
    x3_ = x3 + K3*(z3-H*x3);            % Estimate current state
    %P3_ = (eye(6)-K3*H)*P3*(eye(6) -K3*H)' + K3*R*K3';  % Update the estimate uncertainty
    %x4 = A*x3_;                         % Predict next state
    %P4 = A*P3_*A' + Q;                  % Predict next state

    % Work out the intital collision vertex
    PVx(i) = x3_(1) - (x3_(5)*x3_(2)/x3_(6));
    PVy(i) = x3_(3) - (x3_(5)*x3_(4)/x3_(6));
end

end