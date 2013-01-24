function [Ap q R]=accel2quat(Ar,bias)
    ref=3300;
    sensitivity=300;
    scale_factor=ref/(1023*sensitivity);
    Ap=(Ar-bias)*scale_factor;
    Ap(1)=-Ap(1); Ap(2)=-Ap(2);
    %Ap=Ap/sqrt(sum(Ap.*Ap));
    Ap=normc(Ap);
    theta=atan2(Ap(1),sqrt(Ap(2)*Ap(2) + Ap(3)*Ap(3)));
    si=atan2(Ap(2),Ap(3));
    q=angle2quat(0,theta,si);
    R=quat2dcm(q);
    q=q';
end