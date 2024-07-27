%% 
clc;
clear;
%% ����
syms x1 x2; 
%����һ��Rosenbrock����������������ѻ����㷨���ܵķ�͹������Ҳ����Ϊ�㽶�������ú������ص���ֻ��һ��ȫ����͵㡣
f = 100*(x2-x1^2)^2+(1-x1)^2;  
f = matlabFunction(f); % �����ű��ʽת��Ϊ����������ļ�
%% ��ͼ
% ����������ˮƽ����Ȼ��hold on������������ͼ����ʾ��ÿһ�����������͵������
twod = figure(1); clf; fc = fcontour(f,[-1.5 1.5 -1.5 1.5]); axis equal; hold on
fc.LineWidth = 1;
fc.LineStyle = '--';
fc.LevelList = [1200 1000 800 600 400 200 0];
colorbar;
xlabel('$x_1$', 'Interpreter', 'latex', 'fontsize', 10.5);
ylabel('$x_2$', 'Interpreter', 'latex', 'fontsize', 10.5, 'rotation', pi/2);
title('$f(\mathbf{x}) = 100\times (x_2-x_1^2)^2+(1-x_1)^2$', 'Interpreter', 'latex', 'fontsize', 10.5);
%% ��ʼ����
f_sym = sym(f); %���������ת���ط��ű��ʽ
F_td = matlabFunction(gradient(f_sym));% ��һ�ε���
x0 = [-3/4 1]'; % ��һ����������ʼ��
eps = 1e-5; % ��������Ϊһ���ܽӽ�0����
k=0; % ��k=0, ������¼��������
f_td = F_td(x0(1),x0(2)); %����һ�ε�������ʼ���ֵ
% 3d��ͼ������
threedx1 = [x0(1)];threedx2 = [x0(2)];threedf = [f(x0(1),x0(2))];
%% ��ʼ����
while norm(f_td) > eps   %����һ������Ϊһ�ε�����ֵ���ڿ�������ѭ��
    d= - f_td; % �����½����������½����Ķ��壩
    [alpha_armijo] = armijo(x0,d,f,F_td); %Ƕ��Armijo���ѷ�
    xnew = x0 + alpha_armijo*d; %����ʼ����ͨ��Armijo�õ��Ĳ����������µĵ㡣
     % ������hold on��ˮƽ�����滭���µĵ���������㣩��
    plot([x0(1) xnew(1)],[x0(2) xnew(2)],'ko-')
    refresh
    x0 = xnew; %�õ����㸲����ʼ�㣬ÿ��ѭ���͸�����һ����
    threedx1 = [threedx1,x0(1)];threedx2 = [threedx2,x0(2)];threedf = [threedf,f(x0(1),x0(2))];
    f_td = F_td(x0(1),x0(2)); %�����������һ�ε��ڵ����������ֵ������ѭ���������ж�
    disp(['�������� k = ',num2str(k), ',������xk = (', num2str(x0(1)), num2str(x0(2)), ')']);
    k = k+1; %��������+1
end
x = x0; %���ѭ�������õ������յ�������Ǻ��������ŵ�
result = f(x(1),x(2)); %�����ڸ����ŵ��ֵ���ֲ���Сֵ��
hold off;
%% 3d��ͼ
threed = figure(2);
% ����x��y�ķ�Χ
x = linspace(-1.5, 1.5, 25);
y = linspace(-1.5, 1.5, 25);
% ����x��y������
[X, Y] = meshgrid(x, y);
% ����zֵ������ʹ�õ���һ��ɽ����״�ĺ���
Z  = f(X,Y);
% ���Ƶ�ֵ��ͼ
mesh(X,Y,Z);   % ������άͼ
hold on;
axis square;  % ������֮��ĳ߶����
plot3(threedx1, threedx2, threedf, 'ko-');
xlabel('$x_1$', 'Interpreter', 'latex', 'fontsize', 10.5);
ylabel('$x_2$', 'Interpreter', 'latex', 'fontsize', 10.5, 'rotation', pi/2);
zlabel('$f(\mathbf{x})$', 'Interpreter', 'latex', 'fontsize', 10.5);
title('$f(\mathbf{x}) = 100\times (x_2-x_1^2)^2+(1-x_1)^2$', 'Interpreter', 'latex', 'fontsize', 10.5);