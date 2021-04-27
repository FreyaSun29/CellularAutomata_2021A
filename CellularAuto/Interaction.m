clear,clc
% �����ʼ����
% m��n������100*100��Ԫ��
% p�������
% h�����������
m=100;
n=100;
p1=.7; %AB������B�ĸ���
p2=.15;%AB������A�ĸ���
p3=.15;%AB������յصĸ���
h=122;
% A_rate = 1;%A����������
% B_rate = 2;%B����������
% C_rate = 3;%C����������

% 1.���ɳ�ʼԪ��
% ��ѭ������ÿһ��Ԫ��
for x=1:m
    for y=1:n 
        % ����һ��0~1�����
        r=rand(1);
        % p = 0.3,���������0~0.3֮�俴��Ϊ��ʼ��0.3~1�����յ�
        if r<0.005
            a(x,y)=1; %A
        elseif r>0.005 && r<0.01
            a(x,y)=2; %B
        elseif r>0.01 && r<0.015
            a(x,y)=3; %C    
        else
            a(x,y)=0; %�յ�
        end
    end
end

% 2.���ɳ�ʼԪ��ͼ
for x=1:m 
    for y=1:n
        % ��ȫ��Ѱ��Ϊ1�ģ�A����Ϳ�ϻ�ɫ��2Ϳ����ɫ��3Ϳ�Ϻ�ɫ
        % ���ھ�����˵�Ǵ�(1,1)��ʼ�ģ���ô(1,1)�Ϳ�����ͼ��
        % (1,1),(1,0),(0,1),(0,0)�ķ����������(x-1,y-1)
        % (x-1,y),(x,y-1),(x,y)�ĸ�����ɵķ����������е�
        % (x,y)
        if a(x,y)==1   
            fx=[x-1,x-1,x,x];
            fy=[y-1,y,y,y-1];
            fill(fx,fy,'y')
            hold on
            % ��������һ��ƴ�Ӿ���һ����ŵ�������ھ�ģ����
        elseif a(x,y)==2
            fx=[x-1,x-1,x,x];
            fy=[y-1,y,y,y-1];
            fill(fx,fy,'b')
            hold on
        elseif a(x,y)==3
            fx=[x-1,x-1,x,x];
            fy=[y-1,y,y,y-1];
            fill(fx,fy,'r')
            hold on
        end
    end
end

% 3.���ݹ������122��
for k=1:h
    % ��¼��ʷ�仯������������¼122��Ԫ������
    % ��Ҫ�õ���ά����
    history(:,:,k) = a;
    % (1)�Ȱ�ȫ����Ϳ�ڣ��ٸ�������λ�û���ɫ
    % ��Ȼͼ��ȫ�������ˣ����Ǿ���û�б��ı䣬��˲��ص���
    % ��Ϣȫ����ʧ
%     fx=[0,m,m,0];
%     fy=[0,0,n,n];
%     fill(fx,fy,'k')
%     hold on
    % (2)��ʼ���ݹ������Ԫ������
    % ��cΪ�¾���aΪԭ����
    for x=2:m-1 
        for y=2:n-1
            % ������Χ�м���A,B,C
            A = 0;
            B = 0;
            C = 0;
            %1
            if a(x-1,y-1) == 1
                A = A+1;
            elseif a(x-1,y-1) == 2
                B = B+1;
            elseif a(x-1,y-1) == 3
                C = C+1;
            end
            %2
            if a(x-1,y) == 1
                A = A+1;
            elseif a(x-1,y) == 2
                B = B+1;
            elseif a(x-1,y) == 3
                C = C+1;
            end
            %3
            if a(x-1,y+1) == 1
                A = A+1;
            elseif a(x-1,y+1) == 2
                B = B+1;
            elseif a(x-1,y+1) == 3
                C = C+1;
            end
            %4
            if a(x,y-1) == 1
                A = A+1;
            elseif a(x,y-1) == 2
                B = B+1;
            elseif a(x,y-1) == 3
                C = C+1;
            end
            %5
            if a(x,y+1) == 1
                A = A+1;
            elseif a(x,y+1) == 2
                B = B+1;
            elseif a(x,y+1) == 3
                C = C+1;
            end
            %6
            if a(x+1,y-1) == 1
                A = A+1;
            elseif a(x+1,y-1) == 2
                B = B+1;
            elseif a(x+1,y-1) == 3
                C = C+1;
            end
            %7
            if a(x+1,y) == 1
                A = A+1;
            elseif a(x+1,y) == 2
                B = B+1;
            elseif a(x+1,y) == 3
                C = C+1;
            end
            %8
            if a(x+1,y+1) == 1
                A = A+1;
            elseif a(x+1,y+1) == 2
                B = B+1;
            elseif a(x+1,y+1) == 3
                C = C+1;
            end
            % ��һ���������пյ�������Χ�ھ�����������Χû���ھ�
            % ��������������һ������
            if a(x,y) == 0
                if A==0 && B==0 && C==0
                    c(x,y) = a(x,y);
                elseif A~=0 && B==0 && C==0
                    r = rand(1);
                    if r<1/8
                        c(x,y) = 1;
                    else
                        c(x,y) = 0;
                    end
                elseif A==0 && B~=0 && C==0
                    r = rand(1);
                    if r<3/8
                        c(x,y) = 2;
                    else
                        c(x,y) = 0;
                    end
                elseif A==0 && B==0 && C~=0
                    r = rand(1);
                    if r<6/8
                        c(x,y) = 3;
                    else
                        c(x,y) = 0;
                    end
                elseif A==0 && B~=0 && C~=0
                    r = rand(1);
                    if r<4/8
                        c(x,y) = 3;
                    elseif r<3/8
                        c(x,y) = 2;
                    else
                        c(x,y) = a(x,y);
                    end
                elseif A~=0 && B==0 && C~=0
                    r = rand(1);
                    if r<3/8
                        c(x,y) = 1;
                    elseif r<3/8
                        c(x,y) = 3;
                    else
                        c(x,y) = a(x,y);
                    end
                elseif A~=0 && B~=0 && C==0
                    r = rand(1);
                    if r<2/8
                        c(x,y) = 1;
                    elseif r<4/8
                        c(x,y) = 2;
                    else
                        c(x,y) = a(x,y);
                    end
                elseif A~=0 && B~=0 && C~=0
                    r = rand(1);
                    if r<1/8
                        c(x,y) = 1;
                    elseif r<3/8
                        c(x,y) = 2;
                    elseif r<3/8
                        c(x,y) = 3;
                    else
                        c(x,y) = a(x,y);
                    end
                end
            end
            % �ڶ������򣬶���A��˵����Χ��B���и��ʱ��ı䣬��C��
            % �仯,��A�ޱ仯
            if a(x,y) == 1
                % ������ϣ��Ը�����ʽ�жϸõ�
                % ����A��˵ֻ������B�ſ��ܱ仯
                if B ~= 0
                    r = rand(1);
                    if r<0.4 %0.4���ʱ��B
                        c(x,y) = 2;
                    elseif r>0.4 && r<0.7 %0.3���ʲ���
                        c(x,y) = a(x,y);
                    elseif r>0.7 && r<1 %0.3��յ�
                        c(x,y) = 0;
                    end
                else 
                    c(x,y) = a(x,y);
                end
            end
            % ���������򣬶�B��˵����Χ��A���ܱ��ģ���C���ܱ���
            if a(x,y) == 2
                if A~=0 && C == 0
                    r = rand(1);
                    if r<0.4 %0.4���ʲ���
                        c(x,y) = a(x,y);
                    elseif r>0.4 && r<0.7 %0.3���ʱ��A
                        c(x,y) = 1;
                    elseif r>0.7 && r<1 %0.3��յ�
                        c(x,y) = 0;
                    end
                elseif A==0 && C ~= 0
                    r = rand(1);
                    if r<0.33 %0.33���ʲ���
                        c(x,y) = a(x,y);
                    elseif r>0.33 && r<0.66 %0.33���ʱ��C
                        c(x,y) = 3;
                    elseif r>0.66 && r<1 %0.33��յ�
                        c(x,y) = 0;
                    end
                elseif A~=0 && C ~= 0
                    r = rand(1);
                    if r<0.45 %0.45���ʲ���
                        c(x,y) = a(x,y);
                    elseif r>0.45 && r<0.6 %0.15���ʱ��A
                        c(x,y) = 1;
                    elseif r>0.6 && r<0.75 %0.15���ʱ��C
                        c(x,y) = 3;
                    elseif r>0.66 && r<1 %0.33��յ�
                        c(x,y) = 0;
                    end
                else 
                    c(x,y) = a(x,y);
                end
            end
            % ���������򣬶�C��˵����Χ��A���䣬��B���ܱ���
            if a(x,y) == 3
                if  B~=0
                    r = rand(1);
                    if r<0.33 %0.33���ʱ��B
                        c(x,y) = 2;
                    elseif r>0.33 && r<0.66 %0.33���ʲ���
                        c(x,y) = a(x,y);
                    elseif r>0.66 && r<1 %0.3��յ�
                        c(x,y) = 0;
                    end
                else
                    c(x,y) = a(x,y);
                end
            end
        end
    end
    % û�������ھ��򱣳�ԭ����
    c(1:m,1)=a(1:m,1);c(1:m,n)=a(1:m,n);
    % (3)����Ԫ������ɫ
    for x=1:m 
        for y=1:n
            if c(x,y)==1 
                fx=[x-1,x-1,x,x];
                fy=[y-1,y,y,y-1];
                fill(fx,fy,'y'),hold on
            elseif c(x,y)==2 
                fx=[x-1,x-1,x,x];
                fy=[y-1,y,y,y-1];
                fill(fx,fy,'b'),hold on
            elseif c(x,y)==3 
                fx=[x-1,x-1,x,x];
                fy=[y-1,y,y,y-1];
                fill(fx,fy,'r'),hold on
            end
        end
    end
    % ����һ��
    pause(0.1)
    % �־����Ϊԭ���󣬽�����һ��ѭ��
    a=c;
end

% % 4.����history��������������ͼ
% live = [];
% dead = [];
% for k=1:h
%     live(k) = sum(sum(history(:,:,k) == 1));
%     dead(k) = sum(sum(history(:,:,k) == 0));
% end
% 
% plot(1:h,live,"g",1:h,dead,"k");
% legend("live","dead")