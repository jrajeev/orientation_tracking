%% Load datas
%loading vicon data
load('../data/viconRot1.mat');
tsvicon=ts;
%loading imu data
load('../data/imuRaw1.mat');
tsimu=ts;
%loading cam data
load('../data/cam1.mat');
tscam=ts;

%% Synchronize time stamps
[tsimu_vicon]=synchronize_ts(tsimu,tsvicon);
[tsvicon_cam]=synchronize_ts(tsvicon,tscam);
[tsimu_cam]=synchronize_ts(tsimu,tscam);

%% Convert data to states and observations
bias_accel=[510.79; 501; 512.13];
[accel,R_accel]=accel2rotmat(vals(1:3,:),bias_accel);
bias_gyro=[369.66; 373.63; 375.2];
[gyro,R_gyro]=gyro2rotmat(vals(4:6,:),bias_gyro);

%% Unscented Kalman Filter - 7 states
P=diag(rand(1,6));
Q=diag(rand(1,6));
R=diag(rand(1,6));
[q_fin, R_ukf]=ukf_7state(accel,gyro,P,Q,R);

%% Plot the results & error metric
for i=1+tsimu_vicon:numel(tsimu)
    result_plot(R_ukf(:,:,i),rots(:,:,i-tsimu_vicon),R_accel(:,:,i),R_gyro(:,:,i),i);
end
%% To do
% How to compute mean
% EM algorithm for fixing Qq and R