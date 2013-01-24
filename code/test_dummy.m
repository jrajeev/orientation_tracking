%% Load files
%loading vicon data
load('../data/viconRot1.mat');
tsvicon=ts;
%loading imu data
load('../data/imuRaw1.mat');
tsimu=ts;
%loading cam data
load('../data/cam1.mat');
tscam=ts;

%%
P=diag(rand(1,6));
Q=diag(rand(1,6));
R=diag(rand(1,6));
del_t=1/100; %time step
    %states initialized from first accelerometer reading
    grad_desc_thresh=0.01;
    g=[0 0 0 1];
    R_ukf=zeros(3,3,size(accel,2));
    X=zeros(7,12);
    Y=zeros(7,12);
    Z=zeros(6,12);
    omega=[0;0;0];
    xcap=[1;0;0;0;0;0;0];
    xcap_=[1;0;0;0;0;0;0];
    q_fin=zeros(4,size(accel,2));
    for i=1:10%size(accel,2)
        % Quaternion sigma points
        fprintf('Timestep # %d\n',i);
        L=chol(P+Q);
        L1=L*sqrt(12); %n=6 sqrt(2n)
        L2=-L1;
        W=[L1 L2];
        X(1:4,:)=vect2quat(W(1:3,:));
        X(5:7,:)=W(4:6,:);
        %Process Model
        % Transformation of the sigma points
        norm_factor=sqrt(sum(omega.*omega));
        alpha_del=norm_factor*del_t;
        e_del=omega/norm_factor;
        e_del(isnan(e_del))=0;
        q_del=[cos(alpha_del/2); e_del*sin(alpha_del/2)];

        q_int=quatmultiply(transpose(xcap(1:4)),transpose(X(1:4,:)));
        qi=quatmultiply(q_int,q_del');
        qi=quatnormalize(qi);
        Y(1:4,:)=qi';
        Y(5:7,:)=bsxfun(@plus,omega,X(5:7,:));
        % Gradient descent for a priori mean computation
        q=transpose(xcap(1:4));
        iters=0;
        norm_e=1;
        while (norm_e>grad_desc_thresh &&  iters<50)
            e_quat=quatmultiply(qi,quatinv(q));
            e_vect=quat2vect(e_quat');
            e_mean=mean(e_vect,2);
            e_quat=vect2quat(e_mean);
            e_quat=e_quat';
            q=quatmultiply(e_quat,q);
            norm_e=sqrt(sum(e_mean.*e_mean));
            iters=iters+1;
            %fprintf('iteration # %d\n',iters);
        end
        fprintf('gradient descent algo took  %d  iterations and norm = %d\n',iters,norm_e);
        xcap_(1:4)=q';
        xcap_(5:7)=mean(Y(5:7,:),2);
        
        % Computation of a priori state vector covariance
        Wprime=[e_vect;Y(5:7,:)];
        P_=cov(Wprime');
        
        % Measurement Estimate Covariance
        gprime=quatmultiply(transpose(Y(1:4,:)),quatmultiply(g,quatinv(transpose(Y(1:4,:)))));
        gvect=quat2vect(gprime'); 
        Z(1:3,:)=gvect;
        Z(4:6,:)=Y(5:7,:)
        z_=mean(Z,2);
        z_obs=[accel(:,i);gyro(:,i)];
        v=z_obs-z_
        Pzz=cov(Z')
        Pvv=Pzz + R;

        % Cross Correlation Matrix
        tmp=bsxfun(@minus,Z,z_);
        tmp=tmp';
        Pxz=zeros(6,6);
        for j=1:size(Wprime,2)
            Pxz=Pxz+(Wprime(:,j)*tmp(j,:));
        end
        Pxz=Pxz/6;
        
        %Update equations
        K=Pxz*(Pvv^-1); % Kalman gain update
        tmp=K*v;
        xcap(1:4)=transpose(quatnormalize(quatmultiply(transpose(xcap_(1:4)),transpose(vect2quat(tmp(1:3)))))); % a posteriori estimate update
        xcap(5:7)=xcap_(5:7)+tmp(4:6);
        %P=P_ - K*Pvv*K'; % State covariance update;
        P=P_;
        %Convert xcap to quaternion and then to rotation matrix
        q_fin(:,i)=xcap(1:4);
        omega=xcap(5:7);
        R_ukf(:,:,i)=quat2dcm(transpose(q_fin(:,i)));
    end