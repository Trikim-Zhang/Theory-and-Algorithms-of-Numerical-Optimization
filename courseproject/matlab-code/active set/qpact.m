function [x,lamk,exitflag,output, threedx1, threedx2]=qpact(H,c,Ae,be,Ai,bi,x0)
% 功能: 用有效集方法解一般约束二次规划问题:  
%   min f(x)=0.5*x'*H*x+c'*x,
%   s.t.a'_i*x-b_i=0,(i=1,...,l),
%       a'_i*x-b_i>=0,(i=l+1,...,m)
% 输入: x0是初始点, H, c分别是目标函数二次型矩阵和向量；
%   Ae=(a_1,...,a_l)',  be=(b_1,...,b_l)'; 
%   Ai=(a_{l+1},...,a_m), bi=(b_{l+1},...,b_m)'.
% 输出:  x是最优解， lambda是对应的乘子向量；output是结构
% 变量, 输出极小值f(x), 迭代次数k等信息, exitflag是算法终止类型
    
    % 绘图
    threedx1 = [x0(1)]; threedx2 = [x0(2)];
    
    % 初始化
    %epsilon 是一个很小的数，用于辅助不等式判断，
    %如在进行 >= 的判断时，往往是 x>=b+epsilon，err也是一个类似的量
    epsilon=1.0e-9; err=1.0e-6;
    k=0; x=x0;  %k为迭代次数，x0为初始点
    n=length(x);  kmax=1.0e3;   %kmax最大迭代次数
    ne=length(be); ni=length(bi); lamk=zeros(ne+ni,1);
    index=ones(ni,1);
    for (i=1:ni)
        if(Ai(i,:)*x>bi(i)+epsilon), index(i)=0; end
    end
    %算法主程序
    while (k<=kmax)
        disp(['迭代次数 k = ',num2str(k), ',迭代点x',num2str(k) ,'= (', num2str(x(1)), ',' ,num2str(x(2)), ')']);
        %求解子问题
        Aee=[ ];
        if(ne>0),  Aee=Ae; end
        for(j=1:ni)        %不等式约束的个数
            if(index(j)>0), Aee=[Aee; Ai(j,:)]; end    %将不等式约束和等式约束的稀疏矩阵合并
        end
    %   min f(x)=0.5*x'*H*x+c'*x,
    %    s.t.  a'_i*x-b_i=0,(i=1,...,l，l+1,...,m),
    %    将不等式约束改为等式约束求解子问题。      
        gk=H*x+c;
        [m1,n1] = size(Aee);
        [dk,lamk]=qsubp(H,gk,Aee,zeros(m1,1));       %计算出极小点dk和拉格朗日乘子向量lamk   
        disp('  求解子问题');
        dk
        lamk
        if(norm(dk)<=err)                            %dk为0的时候转入第二步
            disp(['  d',num2str(k),'=0, ','转入第二步']);
            y=0.0;
            if(length(lamk)>ne)                      
                [y,jk]=min(lamk(ne+1:length(lamk))); %确定有效集中lambda的最小元素。 
            end
            if(y>=0)
                exitflag=0;                          %如果每个lambda都大于零，则dk为全局极小点，
                disp([' 每个λ都大于零，','d',num2str(k),'为全局极小点']);
            else                                     %否则减去lamda对应的有效集元素，形成新的有效集
                exitflag=1;
                for(i=1:ni)
                    if(index(i) && (ne+sum(index(1:i)))==jk)   %如果lambda对应的有效集位置为jk，且索引为1，则将索引置0
                        disp(['  更新有效集S',num2str(k),',去掉i',num2str(k),'=',num2str(i),'和x']);
                        xnew = x;
                        index(i)=0; break;                    %确保在之后的计算中，不在计算当前不等式约束。
                    end
                end
            end
%             k=k+1;
%             disp(['迭代次数 k = ',num2str(k), ',迭代点x',num2str(k) ,'= (', num2str(x(1)), ',' ,num2str(x(2)), ')']);
        else                                         %如果dk不等于0，转入第三步
            disp(['  d',num2str(k),'≠0, ','转入第三步，确定步长α']);
            exitflag=1;                              %确定步长alpha
            %求步长
            alpha=1.0; tm=1.0;
            for(i=1:ni)
                if((index(i)==0)&(Ai(i,:)*dk<0))
                    tm1=(bi(i)-Ai(i,:)*x)/(Ai(i,:)*dk);  %alpha的计算见第三步的具体公式
                    if(tm1<tm)
                        tm=tm1; ti=i;
                    end
                end
            end
            disp('  更新有效集和x');
            alpha=min(alpha,tm);                        %选取最小的alpha，一般为边界点对应的alpha
            xnew = x+alpha*dk;                             %确定新的x位置。
            %修正有效集
            if(tm<1), index(ti)=1; end
        end
        if(exitflag==0), break; end
        plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-');
        refresh;
        x = xnew;
        threedx1 = [threedx1;x(1)]; threedx2 = [threedx2;x(2)];
        k=k+1;
    end

    output.fval=0.5*x'*H*x+c'*x;
    output.iter=k;
end