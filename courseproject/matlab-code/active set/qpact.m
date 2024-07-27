function [x,lamk,exitflag,output, threedx1, threedx2]=qpact(H,c,Ae,be,Ai,bi,x0)
% ����: ����Ч��������һ��Լ�����ι滮����:  
%   min f(x)=0.5*x'*H*x+c'*x,
%   s.t.a'_i*x-b_i=0,(i=1,...,l),
%       a'_i*x-b_i>=0,(i=l+1,...,m)
% ����: x0�ǳ�ʼ��, H, c�ֱ���Ŀ�꺯�������;����������
%   Ae=(a_1,...,a_l)',  be=(b_1,...,b_l)'; 
%   Ai=(a_{l+1},...,a_m), bi=(b_{l+1},...,b_m)'.
% ���:  x�����Ž⣬ lambda�Ƕ�Ӧ�ĳ���������output�ǽṹ
% ����, �����Сֵf(x), ��������k����Ϣ, exitflag���㷨��ֹ����
    
    % ��ͼ
    threedx1 = [x0(1)]; threedx2 = [x0(2)];
    
    % ��ʼ��
    %epsilon ��һ����С���������ڸ�������ʽ�жϣ�
    %���ڽ��� >= ���ж�ʱ�������� x>=b+epsilon��errҲ��һ�����Ƶ���
    epsilon=1.0e-9; err=1.0e-6;
    k=0; x=x0;  %kΪ����������x0Ϊ��ʼ��
    n=length(x);  kmax=1.0e3;   %kmax����������
    ne=length(be); ni=length(bi); lamk=zeros(ne+ni,1);
    index=ones(ni,1);
    for (i=1:ni)
        if(Ai(i,:)*x>bi(i)+epsilon), index(i)=0; end
    end
    %�㷨������
    while (k<=kmax)
        disp(['�������� k = ',num2str(k), ',������x',num2str(k) ,'= (', num2str(x(1)), ',' ,num2str(x(2)), ')']);
        %���������
        Aee=[ ];
        if(ne>0),  Aee=Ae; end
        for(j=1:ni)        %����ʽԼ���ĸ���
            if(index(j)>0), Aee=[Aee; Ai(j,:)]; end    %������ʽԼ���͵�ʽԼ����ϡ�����ϲ�
        end
    %   min f(x)=0.5*x'*H*x+c'*x,
    %    s.t.  a'_i*x-b_i=0,(i=1,...,l��l+1,...,m),
    %    ������ʽԼ����Ϊ��ʽԼ����������⡣      
        gk=H*x+c;
        [m1,n1] = size(Aee);
        [dk,lamk]=qsubp(H,gk,Aee,zeros(m1,1));       %�������С��dk���������ճ�������lamk   
        disp('  ���������');
        dk
        lamk
        if(norm(dk)<=err)                            %dkΪ0��ʱ��ת��ڶ���
            disp(['  d',num2str(k),'=0, ','ת��ڶ���']);
            y=0.0;
            if(length(lamk)>ne)                      
                [y,jk]=min(lamk(ne+1:length(lamk))); %ȷ����Ч����lambda����СԪ�ء� 
            end
            if(y>=0)
                exitflag=0;                          %���ÿ��lambda�������㣬��dkΪȫ�ּ�С�㣬
                disp([' ÿ���˶������㣬','d',num2str(k),'Ϊȫ�ּ�С��']);
            else                                     %�����ȥlamda��Ӧ����Ч��Ԫ�أ��γ��µ���Ч��
                exitflag=1;
                for(i=1:ni)
                    if(index(i) && (ne+sum(index(1:i)))==jk)   %���lambda��Ӧ����Ч��λ��Ϊjk��������Ϊ1����������0
                        disp(['  ������Ч��S',num2str(k),',ȥ��i',num2str(k),'=',num2str(i),'��x']);
                        xnew = x;
                        index(i)=0; break;                    %ȷ����֮��ļ����У����ڼ��㵱ǰ����ʽԼ����
                    end
                end
            end
%             k=k+1;
%             disp(['�������� k = ',num2str(k), ',������x',num2str(k) ,'= (', num2str(x(1)), ',' ,num2str(x(2)), ')']);
        else                                         %���dk������0��ת�������
            disp(['  d',num2str(k),'��0, ','ת���������ȷ��������']);
            exitflag=1;                              %ȷ������alpha
            %�󲽳�
            alpha=1.0; tm=1.0;
            for(i=1:ni)
                if((index(i)==0)&(Ai(i,:)*dk<0))
                    tm1=(bi(i)-Ai(i,:)*x)/(Ai(i,:)*dk);  %alpha�ļ�����������ľ��幫ʽ
                    if(tm1<tm)
                        tm=tm1; ti=i;
                    end
                end
            end
            disp('  ������Ч����x');
            alpha=min(alpha,tm);                        %ѡȡ��С��alpha��һ��Ϊ�߽���Ӧ��alpha
            xnew = x+alpha*dk;                             %ȷ���µ�xλ�á�
            %������Ч��
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