function [alpha_armijo] = armijo(x,d,f0,f1)
% x是每次循环中得到的函数的解，d是下降方向
% f0是原目标函数， f1是目标函数的一次导数
alpha = 1; % a>0
sigma = 0.4; % 取值范围(0,1)越大越慢
gamma = 0.5; % 取值范围(0,0.5)越大越快
j = 1; %循环标志
    while (j>0)
        x_new = x+alpha.*d;
        if f0(x_new(1),x_new(2))<=f0(x(1),x(2))+sigma*alpha.*f1(x(1),x(2))'*d %检验条件
            j = 0; % 循环break
            alpha_armijo = alpha; 
        else
            alpha = alpha*gamma; %缩小alpha，进入下一次循环
        end    
    end
end
