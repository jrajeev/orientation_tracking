%% Accelerometer vector pointing downwards
%loading imu data
load('../data/imuRaw1.mat');
tsimu=ts;
load('../data/viconRot1.mat');
tsvicon=ts;
bias_accel=[510.79; 501; 512];
[accel, R_accel]=accel2rotmat(vals(1:3,:),bias_accel);
x=[0;0]; y=[0;0]; z=[0;0];
base_accel=zeros(size(accel));
for i=16:size(accel,2)
    base_accel(:,i)=transpose(R_accel(:,:,i))*accel(:,i);
    x(2)=base_accel(1,i); y(2)=base_accel(2,i); z(2)=base_accel(3,i);
    axis([-2 2 -2 2 -2 2])
    axis equal
    hold on
    plot3(x,y,z,'LineWidth',3);
    hold off
    figure,rotplot(rots(:,:,i-15));
    title (['Iteration ',num2str(i)])
    pause(0.01);
end

%%
clear
clc