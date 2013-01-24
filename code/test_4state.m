%%
clear
clc
%% Load datas
%loading vicon data
%load('viconRot11.mat');
%tsvicon=ts;
%q_vicon=dcm2quat(rots);
%q_vicon=q_vicon';
%loading imu data
load('imuRaw9.mat');
tsimu=ts;
%loading cam data
load('cam9.mat');
tscam=ts;

%% Synchronize time stamps
thresh_ts=0.001;
%[tsvicon_imu,tsimu_vicon]=synchronize_ts(tsvicon,tsimu,thresh_ts);
%[tsvicon_cam,tscam_vicon]=synchronize_ts(tsvicon,tscam,thresh_ts);
[tsimu_cam,tscam_imu]=synchronize_ts(tsimu,tscam,thresh_ts);

%% Convert data to states and observations
bias_accel=[510.79; 501; 501];
[Av,q_accel,R_accel]=accel2quat_vect(vals(1:3,:),bias_accel);
bias_gyro=[369.66; 373.63; 375.2];
[Gp, q_gyro, q_delta]=gyro2vect_vect(vals(4:6,:),bias_gyro);

%% Computation of covariance matrix for initialization
cnt=1;
num=numel(find(tsimu_vicon));
err_accel=zeros(3,num);
err_gyro=zeros(3,num);
for i=1:numel(tsimu)
    if (tsimu_vicon(i)>0)
        v_vicon=quat2vect(q_vicon(:,tsimu_vicon(i)));
        v_accel=quat2vect(q_accel(:,i));
        v_gyro=quat2vect(q_gyro(:,i));
        err_accel(:,cnt)=(v_vicon-v_accel);
        err_gyro(:,cnt)=(v_vicon-v_gyro);
        cnt=cnt+1;
    end
end
err_total=[err_accel;err_gyro];
Qq=cov(err_accel');
data_total=[Av; Gp];
P=cov(Av');
R=cov(Av');
%% Unscented Kalman Filter - 4 states
%initialization 
%P=diag(rand(1,3));
%P=cov(transpose(Av(:,100)));
%P=cov(Av');
%R=diag(rand(1,3));
%Qq=diag(rand(1,3));
%P=diag(rand(1,3));%cov of 100 gyro reading
%R=diag(rand(1,3));
%Qq=diag(rand(1,3));
P=diag([0.5 0.4 0.5]);
Qq=diag([0.7 0.6 0.7]);
R=diag([0.3 0.3 0.3]);
del_t=1/100;
zero_mean=[0; 0; 0];
X=zeros(4,6);
%states initialized from first accelerometer reading
q=[1;0;0;0];
thresh=0.0001;
g=[0 0 0 1];
xcap=q;
%rotplot(quat2dcm(transpose(vect2quat(xcap))));
%title('start position');
%pause;
xcap_=zeros(4,1);
Rukf=zeros(3,3,numel(tsimu));
q_fin=zeros(4,numel(tsimu));
tic
for i=1:numel(tsimu)
    % Quaternion sigma points
    fprintf('Timestep # %d\n',i);
    L=chol(P+Qq);
    L1=L*sqrt(6); %n=3 sqrt(2n)
    L2=-L1;
    W=[L1 L2];
    omega=Gp(:,i);
    norm_factor=sqrt(sum(omega.*omega));
    alpha_del=norm_factor*del_t;
    e_del=omega/norm_factor;
    q_del=[cos(alpha_del/2); e_del*sin(alpha_del/2)];
    
    q_W=vect2quat(W);
    q_p=xcap;
    q_int=quatmultiply(q_p',q_W');
    qi=quatmultiply(q_int,q_del');
    %qi=quatnormalize(quatmultiply(q_p',quatmultiply(q_W',q_del')));
    %qi=quatmultiply(q_p',q_W');
    X=qi';
    
    % Transformation of the sigma points
    % Gradient descent for a priori mean computation
    q=q_p';
    %q=quatnormalize(rand(1,4));
    %rotplot(quat2dcm(q));
    %title('before mean');
    %pause;
    iters=1;
    norm_e=1;
    while (norm_e>thresh &&  iters<50)
        fprintf('iteration # %d\n',iters);
        e_quat=quatmultiply(qi,quatinv(q));
        e_vect=quat2vect(e_quat');
        %norm_e=sum(sqrt(sum(e_vect.*e_vect)));
        e_mean=mean(e_vect,2);
        norm_e=sqrt(sum(e_mean.*e_mean));
        e_quat=vect2quat(e_mean);
        e_quat=e_quat';
        q=quatmultiply(e_quat,q);
        iters=iters+1;
    end
    fprintf('gradient descent algo took  %d  iterations and norm = %d\n',iters-1,norm_e);
    %rotplot(quat2dcm(q));
    %title('after mean');
    %pause;
    q=q';
    xcap_=q;
    
    % Computation of a priori state vector covariance
    Wprime=e_vect;
    P_=cov(Wprime');
    %P_=cov(W');
    % Measurement Estimate Covariance
    %vrot=mvnrnd(0,omegacov,size(omega_p,2));
    %zrot=omega_p+ vrot';
    gprime=quatmultiply(X',quatmultiply(g,quatinv(X')));
    gvect=quat2vect(gprime'); 
    %vacc=mvnrnd(zero_mean,R,size(X,2));
    %zacc=gvect+vacc';
    zacc=gvect;
    Z=zacc;
    z_=mean(Z,2);
    z_obs=Av(:,i);
    %g=[0 z_obs']; %for next timestep
    v=z_obs-z_;
    Pzz=cov(Z');
    Pvv=Pzz + R;
    
    % Cross Correlation Matrix
    tmp=bsxfun(@minus,Z,z_);
    tmp=tmp';
    Pxz=zeros(3,3);
    for j=1:size(W,2)
        Pxz=Pxz+(W(:,j)*tmp(j,:));
    end
    Pxz=Pxz/6;
    %sig=1e-3;
    %Pxz=Pxz+sig*eye(3);
    %Update equations
    K=Pxz*(Pvv^-1); % Kalman gain update
    xcap=transpose(quatmultiply(xcap_',transpose(vect2quat(K*v)))); % a posteriori estimate update
    %P=P_ - K*Pvv*K'; % State covariance update;
    %P=P_;
    %Convert xcap to quaternion and then to rotation matrix
    q_fin(:,i)=xcap;
    Rukf(:,:,i)=quat2dcm(transpose(q_fin(:,i)));
    %result_plot(Rukf(:,:,i),rots(:,:,i+15),R_accel(:,:,i),quat2dcm(transpose(q_gyro(:,i))));
    %title(['Iteration ',num2str(i),' / ',num2str(numel(Av,2))]);
    %rotplot(Rukf(:,:,i));
    %title('final ukf');
    %pause;
end
toc
%% Plot the results & error metric
err=zeros(size(tsimu));
for i=1:numel(tsimu)
    if (tsimu_vicon(i)>0)
        tmp=q_vicon(:,tsimu_vicon(i))-q_fin(:,i);
        err(i)=quatnorm(tmp');
        result_plot(Rukf(:,:,i),rots(:,:,tsimu_vicon(i)),R_accel(:,:,i),quat2dcm(transpose(q_gyro(:,i))));
        title(['Iteration ',num2str(i)]);
        fprintf('Vicon has valid timestamp\n');
    else
        result_plot(Rukf(:,:,i),R_accel(:,:,i),R_accel(:,:,i),quat2dcm(transpose(q_gyro(:,i))));
        title(['Iteration ',num2str(i)]);
        fprintf('No valid Vicon timestamp\n');
    end
    pause (0.001);
end

%% Plot the results for valid vicon data
err=zeros(size(tsimu));
idx=find(tsimu_cam);
for i=1:numel(idx)
    %if (tsimu_vicon(i)>0)
        %tmp=q_vicon(:,tsimu_vicon(i))-q_fin(:,i);
        %err(i)=quatnorm(tmp');
        %result_plot(Rukf(:,:,idx(i)),rots(:,:,tsimu_vicon(idx(i))),R_accel(:,:,idx(i)),quat2dcm(transpose(q_gyro(:,idx(i)))));
        result_plot1(Rukf(:,:,idx(i)),cam(:,:,:,tsimu_cam(idx(i)),R_accel(:,:,idx(i)),quat2dcm(transpose(q_gyro(:,idx(i))))));
        title(['Iteration ',num2str(i),' / ',num2str(numel(idx))]);
    %    fprintf('Vicon has valid timestamp\n');
    %else
    %    result_plot(Rukf(:,:,idx(i)),rots(:,:,tsimu_vicon(idx(i))),R_accel(:,:,idx(i)),quat2dcm(transpose(q_gyro(:,idx(i)))));
    %    title(['Iteration ',num2str(i)]);
    %    fprintf('No valid Vicon timestamp\n');
    %end
    pause (0.001);
end

%% Plot the results for valid vicon data
%err=zeros(size(tsimu));
%idx=find(tsimu_vicon);
for i=1:5000%numel(idx)
    %if (tsimu_vicon(i)>0)
        %tmp=q_vicon(:,tsimu_vicon(i))-q_fin(:,i);
        %err(i)=quatnorm(tmp');
        %result_plot(Rukf(:,:,i),rots(:,:,i+15),R_accel(:,:,i),quat2dcm(transpose(q_gyro(:,i))));
        result_plot2(Rukf(:,:,i),R_accel(:,:,i),quat2dcm(transpose(q_gyro(:,i))));
        %result_plot1(Rukf(:,:,i),cam(:,:,:,i),R_accel(:,:,i),quat2dcm(transpose(q_gyro(:,i))));
        title(['Iteration ',num2str(i),' / ',num2str(5000)]);
    %    fprintf('Vicon has valid timestamp\n');
    %else
    %    result_plot(Rukf(:,:,idx(i)),rots(:,:,tsimu_vicon(idx(i))),R_accel(:,:,idx(i)),quat2dcm(transpose(q_gyro(:,idx(i)))));
    %    title(['Iteration ',num2str(i)]);
    %    fprintf('No valid Vicon timestamp\n');
    %end
    pause (0.001);
end
%% To do
% How to compute mean
% EM algorithm for fixing Qq and R
%%
i=1984;
result_plot(Rukf(:,:,i),rots(:,:,tsimu_vicon(i)),R_accel(:,:,i),quat2dcm(transpose(q_gyro(:,i))));

%% Mosaic Generation
idx=find(tsimu_vicon);
thresh_ang=50*pi/180;
mosaic=zeros(2000,3000,3);
tic
for i=1:numel(idx)
    [angy,angz,angx]=quat2angle(transpose(q_fin(:,idx(i))));
    angy=angy*180/pi;
    angz=-angz;
    angy=-angy;
    %angx=-angx;
    if (abs(angy)<thresh_ang)    
        I=cam(:,:,:,tsimu_cam(idx(i)));
        I=double(imrotate(I,angx*180/pi));
        shiftr=10*tan(angy);
        shiftc=10*tan(angz);
        shiftr_std=1000-floor(size(I,1)/2);
        shiftc_std=1500-floor(size(I,2)/2);
        sr=floor(shiftr_std+shiftr);
        sc=floor(shiftc_std+shiftc);
        wtmatI=0.5*ones(size(I,1),size(I,2),3);
        wtmatI=(I>0).*wtmatI;
        mosaic(sr:sr+size(I,1)-1,sc:sc+size(I,2)-1,:)=uint8(((1-wtmatI).*mosaic(sr:sr+size(I,1)-1,sc:sc+size(I,2)-1,:)+wtmatI.*I(:,:,:)));
    end
end
G=fspecial('gaussian');
mosaic(:,:,1)=conv2(mosaic(:,:,1),G,'same');
mosaic(:,:,2)=conv2(mosaic(:,:,2),G,'same');
mosaic(:,:,3)=conv2(mosaic(:,:,3),G,'same');
img_idx=find(mosaic);
[row,col,ch]=ind2sub(size(mosaic),img_idx);
xxmin=min(row); xxmax=max(row);
yymin=min(col); yymax=max(col);
zzmin=min(ch); zzmax=max(ch);
toc
figure,imshow(uint8(mosaic(xxmin:xxmax,yymin:yymax,:)));
%% Mosaic Generation
idx=find(tsimu_cam);
thresh_ang=30*pi/180;
mosaic=zeros(2000,3000,3);
for i=1500:2500
    [angx,angy,angz]=quat2angle(transpose(q_fin(:,idx(i))));
    angz=-angz;
    angy=-angy;
    angx=-angx;
    angy=angy*180/pi;
    if (abs(angy)<thresh_ang)    
        I=cam(:,:,:,tsimu_cam(idx(i)));
        I=double(imrotate(I,angx*180/pi));
        shiftr=320*tan(angy);
        shiftc=320*tan(angz);
        shiftr_std=1000-floor(size(I,1)/2);
        shiftc_std=1500-floor(size(I,2)/2);
        sr=floor(shiftr_std+shiftr);
        sc=floor(shiftc_std+shiftc);
        mosaic(sr:sr+size(I,1)-1,sc:sc+size(I,2)-1,:)=uint8((mosaic(sr:sr+size(I,1)-1,sc:sc+size(I,2)-1,:)+I(:,:,:))/2);
    end
end