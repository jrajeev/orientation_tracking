function [Gp,q_gyro,q_del]=gyro2vect_vect(Gr,bias)
    ref=3300;
    sensitivity=3.33*180/(pi);
    scale_factor=ref/(1023*sensitivity);
    Gp=bsxfun(@minus,Gr,bias);
    Gp=Gp*scale_factor;
    tmp=Gp;
    Gp=[tmp(2:3,:); tmp(1,:)];
    del_t=1/100;
    norm_factor=sqrt(sum(Gp.*Gp));
    alpha_del=norm_factor*del_t;
    e_del=bsxfun(@rdivide,Gp,norm_factor);
    q_del=[cos(alpha_del/2); bsxfun(@times,e_del,sin(alpha_del/2))];
    q_del=q_del';
    %qp=vect2quat([0; 0; 1]);
    %qp=qp';
    qp=[1 0 0 0];
    q_gyro=zeros(size(Gp,2),4);
    for i=1:size(Gp,2)
        q_gyro(i,:)=quatnormalize(quatmultiply(qp,q_del(i,:)));
        qp=q_gyro(i,:);
    end
    q_gyro=q_gyro';
end