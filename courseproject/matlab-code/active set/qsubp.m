function [x,lambda]=qsubp(H,c,Ae,be)
%求解子问题即，求解一个等式约束的二次规划，
    ginvH=pinv(H);
    [m,n]=size(Ae);
    if (m>0)
        rb = Ae*ginvH*c + be;
        lambda = pinv(Ae*ginvH*Ae')*rb;
        x = ginvH*(Ae'*lambda-c);
    else
        x = -ginvH*c;
        lambda = zeros(m,1);
    end
end