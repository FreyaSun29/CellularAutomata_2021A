% 模型的修正：主要加入了环境的突变
clear,clc
% 定义初始变量
% m和n代表有100*100个元胞
% p代表概率
% h代表迭代次数
m=100;
n=100;
h=61;
% 天气假设只有晴，阴和雨，对应三种湿度分别是0,0.5,1
weather = xlsread("tem.xlsx");
% 归一化
weather = (weather-min(weather))/(max(weather)-min(weather));
% 离散化
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
        weather(l) = 0; % 晴
    elseif r>0.5 && r<0.75
        weather(l) = 0.5; % 阴
    else
        weather(l) = 1; % 雨
    end
end

% 1.生成初始元胞
% 用循环遍历每一个元胞
for x=1:m
    for y=1:n 
        % 生成一个0~1随机数
        r=rand(1);
        % p = 0.3,随机数落在0~0.015之间看作为初始，0.015~1看作空地
        if r<0.005
            a(x,y)=1; %A
        elseif r>0.005 && r<0.01
            a(x,y)=2; %B
        elseif r>0.01 && r<0.015
            a(x,y)=3; %C    
        else
            a(x,y)=0; %空地
        end
    end
end

% 2.生成初始元胞图
for x=1:m 
    for y=1:n
        % 在全局寻找为1的（A），涂上黄色，2涂上蓝色，3涂上红色
        % 对于矩阵来说是从(1,1)开始的，那么(1,1)就可以用图中
        % (1,1),(1,0),(0,1),(0,0)的方格代表，即用(x-1,y-1)
        % (x-1,y),(x,y-1),(x,y)四个点组成的方格代表矩阵中的
        % (x,y)
        if a(x,y)==1   
            fx=[x-1,x-1,x,x];
            fy=[y-1,y,y,y-1];
            fill(fx,fy,'y')
            hold on
            % 这两个在一起拼接就是一个冯诺伊曼的邻居模型了
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
% 保持一秒
pause(0.1)

% 3.根据规则迭代122次
for k=1:h
    disp(k)
    % 记录历史变化，可以用来记录122次元胞矩阵
    % 需要用到三维数组
    history(:,:,k) = a;
    % (1)先把全部都涂黑，再根据生的位置画黄色
    % 虽然图被全部画黑了，但是矩阵还没有被改变，因此不必担心
    % 信息全部丢失
%     fx=[0,m,m,0];
%     fy=[0,0,n,n];
%     fill(fx,fy,'k')
%     hold on
    % (2)开始根据规则更新元胞矩阵
    % 以c为新矩阵，a为原矩阵
    for x=2:m-1 
        for y=2:n-1
            % 计算周围有几个A,B,C
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
            % 第一条规则，若有空地则向周围邻居生长，若周围没有邻居
            % 则不生长，继续下一条规则
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
            % 第二条规则，对于A来说，周围有B会有概率被改变，有C无
            % 变化,有A无变化
            if a(x,y) == 1
                % 计算完毕，以概率形式判断该点
                % 对于A来说只有遇到B才可能变化
                if B ~= 0
                    r = rand(1);
                    if r<0.4 %0.4概率变成B
                        c(x,y) = 2;
                    elseif r>0.4 && r<0.7 %0.3概率不变
                        c(x,y) = a(x,y);
                    elseif r>0.7 && r<1 %0.3变空地
                        c(x,y) = 0;
                    end
                else 
                    c(x,y) = a(x,y);
                end
            end
            % 第三条规则，对B来说，周围有A可能被改，有C可能被改
            if a(x,y) == 2
                if A~=0 && C == 0
                    r = rand(1);
                    if r<0.4 %0.4概率不变
                        c(x,y) = a(x,y);
                    elseif r>0.4 && r<0.7 %0.3概率变成A
                        c(x,y) = 1;
                    elseif r>0.7 && r<1 %0.3变空地
                        c(x,y) = 0;
                    end
                elseif A==0 && C ~= 0
                    r = rand(1);
                    if r<0.33 %0.33概率不变
                        c(x,y) = a(x,y);
                    elseif r>0.33 && r<0.66 %0.33概率变成C
                        c(x,y) = 3;
                    elseif r>0.66 && r<1 %0.33变空地
                        c(x,y) = 0;
                    end
                elseif A~=0 && C ~= 0
                    r = rand(1);
                    if r<0.45 %0.45概率不变
                        c(x,y) = a(x,y);
                    elseif r>0.45 && r<0.6 %0.15概率变成A
                        c(x,y) = 1;
                    elseif r>0.6 && r<0.75 %0.15概率变成C
                        c(x,y) = 3;
                    elseif r>0.66 && r<1 %0.33变空地
                        c(x,y) = 0;
                    end
                else 
                    c(x,y) = a(x,y);
                end
            end
            % 第四条规则，对C来说，周围有A不变，有B可能被改
            if a(x,y) == 3
                if  B~=0
                    r = rand(1);
                    if r<0.33 %0.33概率变成B
                        c(x,y) = 2;
                    elseif r>0.33 && r<0.66 %0.33概率不变
                        c(x,y) = a(x,y);
                    elseif r>0.66 && r<1 %0.3变空地
                        c(x,y) = 0;
                    end
                else
                    c(x,y) = a(x,y);
                end
            end
        end
    end
            
    % 没有完整邻居则保持原矩阵
    c(1:m,1)=a(1:m,1);c(1:m,n)=a(1:m,n);
    % 第五条规则，根据天气的情况来选取一些生长快的死亡
    % 当天气的变化是0.5时
    if k>=2
        if abs(weather(k)-weather(k-1))==1 % 差别大，死的多一点
            % 找到B和C的位置
            B_pos = (c==2);
            C_pos = (c==3);
            % 要选择一部分去死，C死的最多
            for x=1:m
                for y=1:n
                    if B_pos(x,y)==1
                        r = rand(1);
                        if r<0.8 %留下0.8的活着
                            B_pos(x,y) = 0;
                        end
                    end
                    if C_pos(x,y)==1
                        r = rand(1);
                        if r<0.6 %留下0.6的活着
                            C_pos(x,y) = 0;
                        end
                    end
                end
            end
            c(B_pos) = 0;
            c(C_pos) = 0;
        elseif abs(weather(k)-weather(k-1))==0.5
            % 找到B和C的位置
            B_pos = (c==2);
            C_pos = (c==3);
            % 要选择一部分去死，C死的最多
            for x=1:m
                for y=1:n
                    if B_pos(x,y)==1
                        r = rand(1);
                        if r<0.9 %留下0.9的活着
                            B_pos(x,y) = 0;
                        end
                    end
                    if C_pos(x,y)==1
                        r = rand(1);
                        if r<0.8 %留下0.8的活着
                            C_pos(x,y) = 0;
                        end
                    end
                end
            end
            c(B_pos) = 0;
            c(C_pos) = 0;
        end
    end
    % (3)根据元胞矩阵画色
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
    % 保持一秒
    pause(0.1)
    % 现矩阵成为原矩阵，进入下一次循环
    a=c;
end

% 4.根据history画三种真菌的趋势图
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
legend("A的数量","B的数量","C的数量")

% 5.根据history求出整块木头被分解的程度
breakdown = ones(m,n); %被分解情况

A_humidity_resistance = -1;%A的生长速率，进入函数换算成死亡率
B_humidity_resistance = 0;%B的生长速率
C_humidity_resistance = 1;%C的生长速率

A_rate = 1;%A的生长速率
B_rate = 2;%B的生长速率
C_rate = 3;%C的生长速率

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

% 去掉首行末行，首列末列
breakdown(:,1) = [];
breakdown(:,end) = [];
breakdown(1,:) = [];
breakdown(end,:) = [];

figure
% 有生物多样性的平均分解程度图
diverse_average_breakdown = (1-sum(sum(breakdown))/size(breakdown,1)/size(breakdown,2))*100
% 根据48x48的breakdown画图
x = 1:size(breakdown,1);
y = 1:size(breakdown,2);
surf(x,y,breakdown);
axis([1,size(breakdown,1),1,size(breakdown,2),0,2])
shading interp;

figure
% 仅有一种生物的多样性图
single_breakdown = ones(size(breakdown,1),size(breakdown,2));
single_breakdown = single_breakdown - 0.005*h;
single_diverse_average_breakdown = (1-sum(sum(single_breakdown))/size(single_breakdown,1)/size(single_breakdown,2))*100
surf(x,y,single_breakdown);
% axis([1,size(single_breakdown,1),1,size(single_breakdown,2),0,2])
shading interp;