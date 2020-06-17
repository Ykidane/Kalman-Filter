%Kalman filter implementation in 2D
close all;
%% parameters to be defined initialy
tic
N = 50;             
t_step = 1;         
%Constant and change in acceleration 
u =[1;1];                     
% u = ones(2,N);               
% %changing
%   acc_change = 25;
%   u(1:acc_change:end) = u(1:acc_change:end)*-1;


x0 =0;              % initial position in x direction
v0 =0;            % initial velocity in x direction
y0=0;            %intial position in y direction
vy0=0;            % intial velocity in y direction

% uncertainties
sigma_x = 1;            % [m] initial process variation standard deviation in position 
sigma_v = 1.3;          % [m/s] initial process variation standard deviation in velocity
sigma_y = 0.8;            % [m] initial process variation standard deviation in position 
sigma_vy= 1.4;          % [m/s] initial process variation standard deviation in velocity
meas_errorx =0.5;      % [m] measurement standard variation in position
meas_errorv = 1.3;      % [m/s] measurmenet standard variation in velocity
meas_errory =1.1;      % [m] measurement standard variation in position
meas_errorvy = 0.7;      % [m/s] measurmenet standard variation in velocity

%  wx=2000; wy=2000; wv=1000; wvy=1000;
%  Q =[wx 0 0 0; 0 wy 0 0; 0 0 wv 0; 0 0 0 wvy]; 
Q=eye(4);
 
t=0:t_step:t_step*(N-1);  % Time vector
S=zeros(4,N);                 % position and velocity in x and y direction. 
S(:,1) = [x0,y0, v0, vy0];



%% Generate measurements

% initial measurements
posx_0 = x0;
velx_0 = v0;
posy_0=y0;
vely_0=vy0;

S_meas = zeros(4,N);   
S_ideal = zeros(4,N);   
   
for j = 2:N
    % position in x direction
    S_ideal(1,j) = posx_0 + v0*t(j) + 0.5*u(1)*t(j)^2;    
    S_meas(1,j) = normrnd(S_ideal(1,j),meas_errorx);  
    
    % velocity in x direction
    S_ideal(3,j) = velx_0 + u(2)*t(j);                    
    S_meas(3,j) = normrnd(S_ideal(3,j),meas_errorv);  
    
        %Position in y direction
    S_ideal(2,j) = posy_0 + vy0*t(j) + 0.5*u(1)*t(j)^2;    
    S_meas(2,j) = normrnd(S_ideal(2,j),meas_errory);  
    
        %Velocity in y direction.
    S_ideal(4,j) = vely_0 + u(2)*t(j);                    
    S_meas(4,j) = normrnd(S_ideal(4,j),meas_errorvy);  
   
end

%% define matrices

% motion
A = [1 0 t_step 0;0 1 0 t_step;0 0 1 0;0 0 0 1];
    % change in position by acceleration  
B = [0.5*t_step^2 0;0  0.5*t_step^2 ;t_step 0;0 t_step];
C = eye(4);
H = eye(4);

% uncertainties
P = [sigma_x^2 0 0 0;0 sigma_y^2 0 0;        
    0 0 sigma_v^2 0;0 0 0 sigma_vy^2];       % off-diagonal elements 0 because  variables are independent
                                                                                                 %observations errors
               
R = [meas_errorx 0 0 0 ;0 meas_errory 0 0;
     0 0 meas_errorv 0;0 0 0 meas_errorvy];             

    
%% Kalman filter

for k=2:N
    % Predicted values
    S_pred = A*S(:,k-1) + B*u ;
  
          
    % Process covariance matrix
    P = A*P*A' + Q;
    P = eye(4).*diag(P);  % omitting covariance terms. 
    
    % Kalman Gain:
    K = P*H'/(H*P*H' + R);
    
    % Transforming measurement values to be compliant with syntax
    S_meas(:,k) = C*S_meas(:,k);
    
    % Current state
    S(:,k) = S_pred + K*((S_meas(:,k)-H*S_pred)); 
    
    % Updating the process covariance matrix
    P = (eye(4)-K*H)*P;
     
end

figure
subplot(2,1,1);
plot(S(1,:), S(2,:), S_meas(1,:),S_meas(2,:),'o',S_ideal(1,:), S_ideal(2,:))
legend({'Kalman filtered','Measured values','Theoretical values'},'Location','northwest')
title(['Position in x direction with sigmx=' num2str(sigma_x),  ' measerrorx=' num2str(meas_errorx), ' and Position in y direction with sigmay=' num2str(sigma_y), ' measerrory=' num2str(meas_errory)])
xlabel('Position in x direction [m]')
ylabel('Position in y direction [m]')


subplot(2,1,2);
plot(S(3,:), S(4,:), S_meas(3,:),S_meas(4,:),'o',S_ideal(3,:), S_ideal(4,:))
legend({'Kalman filtered','Measured values', 'Theoretical values'},'Location','northwest')
title(['Velocity in x direction with sigmv=' num2str(sigma_v),  ' measerrorv=' num2str(meas_errorv), ' and velocity in y direction with sigmavy=' num2str(sigma_vy), ' measerrorvy=' num2str(meas_errorvy)])
xlabel('Velocity in x direction [m/s]')
ylabel('Velocity in y direction [m/s]')
toc

%plot(t,P(2,:))
%  speed=zeros(1,N);
%  for ii=1:N
%      speed(ii)=sqrt(S(3, ii)^2 + S(4,ii)^2); 
%  end
% 
% 
%  figure
% % speed
%  plot(t,speed)
