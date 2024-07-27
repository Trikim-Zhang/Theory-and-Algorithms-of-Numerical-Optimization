%% 
clc;
clear;
%% 函数
syms x1 x2; 
%这是一个Rosenbrock函数，用来测试最佳化演算法性能的非凸函数，也被称为香蕉函数。该函数的特点是只有一个全局最低点。
f = 100*(x2-x1^2)^2+(1-x1)^2;  
f = matlabFunction(f); % 将符号表达式转换为函数句柄或文件
%% 绘图
% 画出函数的水平集，然后hold on，待会会在这个图上显示出每一个迭代步长和迭代结果
twod = figure(1); clf; fc = fcontour(f,[-1.5 1.5 -1.5 1.5]); axis equal; hold on
fc.LineWidth = 1;
fc.LineStyle = '--';
fc.LevelList = [1200 1000 800 600 400 200 0];
colorbar;
xlabel('$x_1$', 'Interpreter', 'latex', 'fontsize', 10.5);
ylabel('$x_2$', 'Interpreter', 'latex', 'fontsize', 10.5, 'rotation', pi/2);
title('$f(\mathbf{x}) = 100\times (x_2-x_1^2)^2+(1-x_1)^2$', 'Interpreter', 'latex', 'fontsize', 10.5);
%% 初始设置
f_sym = sym(f); %将函数句柄转换回符号表达式
F_td = matlabFunction(gradient(f_sym));% 求一次导数
x0 = [-3/4 1]'; % 给一个收敛的起始点
eps = 1e-5; % 设可容误差为一个很接近0的数
k=0; % 设k=0, 用来记录迭代次数
f_td = F_td(x0(1),x0(2)); %计算一次导数在起始点的值
% 3d绘图的数据
threedx1 = [x0(1)];threedx2 = [x0(2)];threedf = [f(x0(1),x0(2))];
%% 开始迭代
while norm(f_td) > eps   %设置一个条件为一次导数的值大于可容误差的循环
    d= - f_td; % 计算下降方向（最速下降法的定义）
    [alpha_armijo] = armijo(x0,d,f,F_td); %嵌套Armijo线搜法
    xnew = x0 + alpha_armijo*d; %在起始点上通过Armijo得到的步长迭代出新的点。
     % 在上面hold on的水平集里面画出新的迭代结果（点）。
    plot([x0(1) xnew(1)],[x0(2) xnew(2)],'ko-')
    refresh
    x0 = xnew; %用迭代点覆盖起始点，每次循环就覆盖上一个点
    threedx1 = [threedx1,x0(1)];threedx2 = [threedx2,x0(2)];threedf = [threedf,f(x0(1),x0(2))];
    f_td = F_td(x0(1),x0(2)); %计算出函数的一次导在迭代点上面的值，用于循环条件的判断
    disp(['迭代次数 k = ',num2str(k), ',迭代点xk = (', num2str(x0(1)), num2str(x0(2)), ')']);
    k = k+1; %迭代次数+1
end
x = x0; %最后循环结束得到的最终迭代点就是函数的最优点
result = f(x(1),x(2)); %函数在该最优点的值（局部最小值）
hold off;
%% 3d绘图
threed = figure(2);
% 定义x和y的范围
x = linspace(-1.5, 1.5, 25);
y = linspace(-1.5, 1.5, 25);
% 创建x和y的网格
[X, Y] = meshgrid(x, y);
% 定义z值，这里使用的是一个山脉形状的函数
Z  = f(X,Y);
% 绘制等值线图
mesh(X,Y,Z);   % 绘制三维图
hold on;
axis square;  % 坐标轴之间的尺度相等
plot3(threedx1, threedx2, threedf, 'ko-');
xlabel('$x_1$', 'Interpreter', 'latex', 'fontsize', 10.5);
ylabel('$x_2$', 'Interpreter', 'latex', 'fontsize', 10.5, 'rotation', pi/2);
zlabel('$f(\mathbf{x})$', 'Interpreter', 'latex', 'fontsize', 10.5);
title('$f(\mathbf{x}) = 100\times (x_2-x_1^2)^2+(1-x_1)^2$', 'Interpreter', 'latex', 'fontsize', 10.5);