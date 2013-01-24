function [gyro,R_gyro]=gyro2rotmat(Gr,bias)
    ref=3300;
    sensitivity=3.33*180/(pi);
    scale_factor=ref/(1023*sensitivity);
    Gp=bsxfun(@minus,Gr,bias);
    Gp=Gp*scale_factor;
    gyro=[Gp(2:3,:); -Gp(1,:)];
    del_t=1/100;
    norm_factor=sqrt(sum(gyro.*gyro));
    alpha_del=norm_factor*del_t;
    e_del=bsxfun(@rdivide,gyro,norm_factor);
    q_del=[cos(alpha_del/2); bsxfun(@times,e_del,sin(alpha_del/2))];
    q_del=q_del';
    qp=[1 0 0 0];
    q_gyro=zeros(size(Gp,2),4);
    for i=1:size(Gp,2)
        q_gyro(i,:)=quatnormalize(quatmultiply(qp,q_del(i,:)));
        qp=q_gyro(i,:);
    end
    R_gyro=quat2dcm(q_gyro);
end