% ģ�͵���������Ҫ�����˻�����ͻ��
clear,clc
% �����ʼ����
% m��n������100*100��Ԫ��
% p�������
% h�����������
m=100;
n=100;
h=61;
% ��������ֻ���磬�����꣬��Ӧ����ʪ�ȷֱ���0,0.5,1
weather = xlsread("tem.xlsx");
% ��һ��
weather = (weather-min(weather))/(max(weather)-min(weather));
% ��ɢ��
for i=1:h
    if weather(i)<0.33
        weather(i) = 0;
    elseif weather(i)>0.33 && weather(i)<0.66
        weather(i) = 0.5;
    else 
        weather(i) = 1;
    end
end
for l=1:h
    r = rand(1);
    if r<0.5
        weather(l) = 0; % ��
    elseif r>0.5 && r<0.75
        weather(l) = 0.5; % ��
    else
        weather(l) = 1; % ��
    end
end

% 1.���ɳ�ʼԪ��
% ��ѭ������ÿһ��Ԫ��
for x=1:m
    for y=1:n 
        % ����һ��0~1�����
        r=rand(1);
        % p = 0.3,���������0~0.015֮�俴��Ϊ��ʼ��0.015~1�����յ�
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
% ����һ��
pause(0.1)

% 3.���ݹ������122��
for k=1:h
    disp(k)
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
    % ���������򣬸��������������ѡȡһЩ�����������
    % �������ı仯��0.5ʱ
    if k>=2
        if abs(weather(k)-weather(k-1))==1 % �������Ķ�һ��
            % �ҵ�B��C��λ��
            B_pos = (c==2);
            C_pos = (c==3);
            % Ҫѡ��һ����ȥ����C�������
            for x=1:m
                for y=1:n
                    if B_pos(x,y)==1
                        r = rand(1);
                        if r<0.8 %����0.8�Ļ���
                            B_pos(x,y) = 0;
                        end
                    end
                    if C_pos(x,y)==1
                        r = rand(1);
                        if r<0.6 %����0.6�Ļ���
                            C_pos(x,y) = 0;
                        end
                    end
                end
            end
            c(B_pos) = 0;
            c(C_pos) = 0;
        elseif abs(weather(k)-weather(k-1))==0.5
            % �ҵ�B��C��λ��
            B_pos = (c==2);
            C_pos = (c==3);
            % Ҫѡ��һ����ȥ����C�������
            for x=1:m
                for y=1:n
                    if B_pos(x,y)==1
                        r = rand(1);
                        if r<0.9 %����0.9�Ļ���
                            B_pos(x,y) = 0;
                        end
                    end
                    if C_pos(x,y)==1
                        r = rand(1);
                        if r<0.8 %����0.8�Ļ���
                            C_pos(x,y) = 0;
                        end
                    end
                end
            end
            c(B_pos) = 0;
            c(C_pos) = 0;
        end
    end
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

% 4.����history���������������ͼ
A_num = [];
B_num = [];
C_num = [];
for k=1:h
    A_num(k) = sum(sum(history(:,:,k) == 1));
    B_num(k) = sum(sum(history(:,:,k) == 2));
    C_num(k) = sum(sum(history(:,:,k) == 3));
end

figure
plot(1:h,A_num,"y",1:h,B_num,"b",1:h,C_num,'r');
legend("A������","B������","C������")

% 5.����history�������ľͷ���ֽ�ĳ̶�
breakdown = ones(m,n); %���ֽ����

A_humidity_resistance = -1;%A���������ʣ����뺯�������������
B_humidity_resistance = 0;%B����������
C_humidity_resistance = 1;%C����������

A_rate = 1;%A����������
B_rate = 2;%B����������
C_rate = 3;%C����������

[t1,x1]=ode45(@(t1,x1)average_breakdown(t1,x1,A_rate,A_humidity_resistance),[1 122],[1 0]);
wood=x1(:,1); fungi=x1(:,2); 
A_average_breakdown = wood(end)/122;
[t2,x2]=ode45(@(t2,x2)average_breakdown(t2,x2,B_rate,B_humidity_resistance),[1 122],[1 0]);
wood=x2(:,1); fungi=x2(:,2); 
B_average_breakdown = wood(end)/122;
[t3,x3]=ode45(@(t3,x3)average_breakdown(t3,x3,C_rate,C_humidity_resistance),[1 122],[1 0]);
wood=x2(:,1); fungi=x2(:,2); 
C_average_breakdown = wood(end)/122;

for k=1:h
    A_logic = (history(:,:,k)==1);
    B_logic = (history(:,:,k)==2);
    C_logic = (history(:,:,k)==3);
    breakdown = breakdown - A_average_breakdown.*A_logic;
    breakdown = breakdown - B_average_breakdown.*B_logic;
    breakdown = breakdown - C_average_breakdown.*C_logic;
end

% ȥ������ĩ�У�����ĩ��
breakdown(:,1) = [];
breakdown(:,end) = [];
breakdown(1,:) = [];
breakdown(end,:) = [];

figure
% ����������Ե�ƽ���ֽ�̶�ͼ
diverse_average_breakdown = (1-sum(sum(breakdown))/size(breakdown,1)/size(breakdown,2))*100
% ����48x48��breakdown��ͼ
x = 1:size(breakdown,1);
y = 1:size(breakdown,2);
surf(x,y,breakdown);
axis([1,size(breakdown,1),1,size(breakdown,2),0,2])
shading interp;

figure
% ����һ������Ķ�����ͼ
single_breakdown = ones(size(breakdown,1),size(breakdown,2));
single_breakdown = single_breakdown - 0.005*h;
single_diverse_average_breakdown = (1-sum(sum(single_breakdown))/size(single_breakdown,1)/size(single_breakdown,2))*100
surf(x,y,single_breakdown);
% axis([1,size(single_breakdown,1),1,size(single_breakdown,2),0,2])
shading interp;