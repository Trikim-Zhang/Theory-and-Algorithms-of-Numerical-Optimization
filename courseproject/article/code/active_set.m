%% 
clc;
clear;
%% 问题
%输入:  x0是初始点, H, c分别是目标函数二次型矩阵和向量；
%   Ae=(a_1,...,a_l)',  be=(b_1,...,b_l)'; 
%   Ai=(a_{l+1},...,a_m), bi=(b_{l+1},...,b_m)'.
%   min f(x)=0.5*x'*H*x+c'*x,
%   s.t. a'_i*x-b_i=0,(i=1,...,l),
%        a'_i*x-b_i>=0,(i=l+1,...,m)
H=[2 -1; -1 4];
c=[-1 -10]';
Ae=[ ]; be=[ ];
Ai=[-3 -2; 1 0; 0 1];
bi=[-6 0 0]';
x0=[0 0]';
%% 绘图
syms x [2,1]; 
f = 0.5*x'*H*x+c'*x;
f = matlabFunction(f); % 将符号表达式转换为函数句柄或文件
twod = figure(1); clf; fc = fcontour(f,[-3 3 -3 3]); axis equal; hold on
fc.LineWidth = 1;
fc.LineStyle = '--';
colorbar;
xlabel('$x_1$', 'Interpreter', 'latex', 'fontsize', 10.5);
ylabel('$x_2$', 'Interpreter', 'latex', 'fontsize', 10.5, 'rotation', pi/2);
title('$f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^{\mathrm{T}}\mathbf{H}\mathbf{x}+\mathbf{c}^{\mathrm{T}}\mathbf{x}$',...
    'Interpreter', 'latex', 'fontsize', 10.5);
%% 求解
[x,lambda,exitflag,output,threedx1, threedx2]=qpact(H,c,Ae,be,Ai,bi,x0)
hold off;
%% 3d绘图
threed = figure(2);
% 定义x和y的范围
x = linspace(-3, 3, 25);
y = linspace(-3, 3, 25);
% 创建x和y的网格
[X, Y] = meshgrid(x, y);
% 定义z值，这里使用的是一个山脉形状的函数
Z  = f(X,Y);
threedf = f(threedx1, threedx2);
% 绘制等值线图
mesh(X,Y,Z);   % 绘制三维图
hold on;
axis square;  % 坐标轴之间的尺度相等
plot3(threedx1, threedx2, threedf, 'ko-');
xlabel('$x_1$', 'Interpreter', 'latex', 'fontsize', 10.5);
ylabel('$x_2$', 'Interpreter', 'latex', 'fontsize', 10.5, 'rotation', pi/2);
zlabel('$f(\mathbf{x})$', 'Interpreter', 'latex', 'fontsize', 10.5);
title('$f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^{\mathrm{T}}\mathbf{H}\mathbf{x}+\mathbf{c}^{\mathrm{T}}\mathbf{x}$',...
    'Interpreter', 'latex', 'fontsize', 10.5);
