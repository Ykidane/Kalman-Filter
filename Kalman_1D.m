%Kalman filter implementation in 2D
tic
N = 50;             
t_step = 1;         

% Constant accelereation and change in acceleration
u=1;
% u = ones(1,N);              % m/s^2 - acceleration 
% acc_change = 25;
% u(1:acc_change:end) = u(1:acc_change:end)*-1;

x0 = 0;              % initial position
v0 = 0;              % initial velocity

% uncertainties
sigma_x = 1;            % [m] initial process variation standard deviation in position 
sigma_v = 1.2;          % [m/s] initial process variation standard deviation in velocity
meas_errorx = 2.1;      % [m] measurement standard variation in position
meas_errorv = 1.5;      % [m/s] measurmenet standard variation in velocity

w_x=1.1;
w_v=1.1;


Q =[w_x  0;
       0 w_v];                   % Process noise covariance matrix

% preallocating etc. 
t=0:t_step:t_step*(N-1);           % time vector
X=zeros(2,N);                         % position and velocity in x. 
                                               % first row: pos. second row: vel. 
X(:,1) = [x0, v0];

%% Generate measurements

% % initial measurements
pos_0 = x0; 
vel_0 = v0;

X_meas = zeros(2,N);    % 1st row = position, 2nd row = velocity
X_ideal = zeros(2,N);   % 1st row = position, 2nd row = velocity

for j = 2:N
    % position    
    X_ideal(1,j) = X_ideal(1,j-1) + X_ideal(2,j-1)*t_step + 0.5*u*t_step^2; % theoretical value
    X_meas(1,j) = normrnd(X_ideal(1,j),meas_errorx);  % measurement / adding noise
    
    % velocity
    X_ideal(2,j) = vel_0 + u*t(j);                    % theoretical values / equation of motion
    X_meas(2,j) = normrnd(X_ideal(2,j),meas_errorv);  % measurement / adding noise
end

% Define matrices

A = [1 t_step;
         0   1   ];          
B = [0.5*t_step^2;
        t_step    ];      
C = eye(2);              % 'maintaining' matrix for measurements. [1 0] -> only measuring position
H = eye(2);

% uncertainties
% process errors
P = [sigma_x^2 0;...        % initial state matrix covariance
    0 sigma_v^2];           % off-diagonal elements 0 bec. variables are independent
% observations errors
R = [meas_errorx 0 ;...
    0 meas_errorv]; % measurement covariance matrix 


% Kalman filter

for k=2:N
    % Predicted values
    X_pred = A*X(:,k-1) + B*u ; 
        
    % Process covariance matrix
    P = A*P*A' + Q;
    P = eye(2).*diag(P);  % omitting covariance terms.  

    % Kalman Gain:
    K = P*H'/(H*P*H' + R);
  
    X_meas(:,k) = C*X_meas(:,k) ; % Transforming measurement values to be compliant with syntax
    In=X_meas(:,k)-H*X_pred;
    
    % Current state
    X(:,k) = X_pred + K*(In); 
    
    % Updating the process covariance matrix
    P = (eye(2)-K*H)*P;
    disp(P);
end

figure
%position
ax1 = subplot(2,1,1);
plot(t,X(1,:),t,X_meas(1,:),'o',t,X_ideal(1,:))
title(['Position with sigmx=' num2str(sigma_x),  ' measerrorx=' num2str(meas_errorx), ' wx=' num2str(w_x),  ' and wv=' num2str(w_v)])
ylabel('Position [m]')
xlabel('Time [s]')
legend({'Kalman filtered','Measured values','Theoretical values'},'Location','northwest')


% velocity
ax2 = subplot(2,1,2);
plot(t,X(2,:),t,X_meas(2,:),'o',t,X_ideal(2,:))
title(['Velocity with sigmav=' num2str(sigma_v), ' measerrorv=' num2str(meas_errorv), ' wx=' num2str(w_x),  ' and wv=' num2str(w_v)])
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend({'Kalman filtered','Measured values','Theoretical values'},'Location','northwest')

toc


