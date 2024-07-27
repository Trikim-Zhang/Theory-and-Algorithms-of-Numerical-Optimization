function [alpha_armijo] = armijo(x,d,f0,f1)
% x��ÿ��ѭ���еõ��ĺ����Ľ⣬d���½�����
% f0��ԭĿ�꺯���� f1��Ŀ�꺯����һ�ε���
alpha = 1; % a>0
sigma = 0.4; % ȡֵ��Χ(0,0.5)Խ��Խ��
gamma = 0.5; % ȡֵ��Χ(0,1)Խ��Խ��
j = 1; %ѭ����־
    while (j>0)
        x_new = x+alpha.*d;
        if f0(x_new(1),x_new(2))<=f0(x(1),x(2))+sigma*alpha.*f1(x(1),x(2))'*d %��������
            j = 0; % ѭ��break
            alpha_armijo = alpha; 
        else
            alpha = alpha*gamma; %��Сalpha��������һ��ѭ��
        end    
    end
end
