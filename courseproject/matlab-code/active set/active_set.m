%% 
clc;
clear;
%% ����
%����:  x0�ǳ�ʼ��, H, c�ֱ���Ŀ�꺯�������;����������
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
%% ��ͼ
syms x [2,1]; 
f = 0.5*x'*H*x+c'*x;
f = matlabFunction(f); % �����ű��ʽת��Ϊ����������ļ�
twod = figure(1); clf; fc = fcontour(f,[-3 3 -3 3]); axis equal; hold on

% ����x=0
x_zero = linspace(-3, 3, 100); % x��ķ�Χ
y_zero = zeros(size(x_zero)); % y���Ӧ��ֵ
plot(x_zero, y_zero, 'k--'); % ����x=0��ֱ��
hold on;
% ����y=0
y_zero = linspace(-3, 3, 100); % x��ķ�Χ
x_zero = zeros(size(y_zero)); % x���Ӧ��ֵ
plot(x_zero, y_zero, 'k--'); % ����y=0��ֱ��
hold on;
% ����ax + by = 0
a = Ai(1,1);b = Ai(1,2);m = bi(1);
y_values = linspace(-3, 3, 100); % x��ķ�Χ
x_values = (m-b*y_values)/a; % ��Ӧ��yֵ
plot(x_values, y_values, 'k--'); % ����ֱ��ax + by = 0

fc.LineWidth = 1;
fc.LineStyle = '--';
colorbar;
xlabel('$x_1$', 'Interpreter', 'latex', 'fontsize', 10.5);
ylabel('$x_2$', 'Interpreter', 'latex', 'fontsize', 10.5, 'rotation', pi/2);
title('$f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^{\mathrm{T}}\mathbf{H}\mathbf{x}+\mathbf{c}^{\mathrm{T}}\mathbf{x}$',...
    'Interpreter', 'latex', 'fontsize', 10.5);
%% ���
[x,lambda,exitflag,output,threedx1, threedx2]=qpact(H,c,Ae,be,Ai,bi,x0)
hold off;
%% 3d��ͼ
threed = figure(2);
% ����x��y�ķ�Χ
x = linspace(-3, 3, 25);
y = linspace(-3, 3, 25);
% ����x��y������
[X, Y] = meshgrid(x, y);
% ����zֵ������ʹ�õ���һ��ɽ����״�ĺ���
Z  = f(X,Y);
threedf = f(threedx1, threedx2);
% ���Ƶ�ֵ��ͼ
mesh(X,Y,Z);   % ������άͼ
hold on;
axis square;  % ������֮��ĳ߶����
plot3(threedx1, threedx2, threedf, 'ko-');
xlabel('$x_1$', 'Interpreter', 'latex', 'fontsize', 10.5);
ylabel('$x_2$', 'Interpreter', 'latex', 'fontsize', 10.5, 'rotation', pi/2);
zlabel('$f(\mathbf{x})$', 'Interpreter', 'latex', 'fontsize', 10.5);
title('$f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^{\mathrm{T}}\mathbf{H}\mathbf{x}+\mathbf{c}^{\mathrm{T}}\mathbf{x}$',...
    'Interpreter', 'latex', 'fontsize', 10.5);
