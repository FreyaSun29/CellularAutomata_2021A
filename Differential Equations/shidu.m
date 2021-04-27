% 不同湿度影响r2，由于是不耐湿的，因此湿度越大，r2越小
% 当其他条件（延申）不变时，只湿度不同的情况下
list = [];
k=1;
for i = 0.5:0.01:0.7
    r1=0.00001;
    r2=i; % 湿度影响r2
    yanshen=1.3;
    lambda1=0.00005;
    lambda2=0.5;
    [t1,x1]=ode45(@(t1,x1)step1(t1,x1,r1,r2,yanshen,lambda1,lambda2),[1 122],[1 0]);
    wood=x1(:,1); fungi=x1(:,2); 
    list(k) = wood(end);
    k = k+1;
end
r2s = 0.5:0.01:0.7;
naishidu = -10*r2s + 6;
plot(naishidu,log(100*(1-list)),"r")
axis([-1.1,1.1,0.5,3.5]);

% 横坐标选为-1到1，代表耐湿度，当耐湿度越接近于1，那么真菌的耐湿程度
% 就越大，那么它死亡的越少，对应的死亡率r2就越小，即-1对应0.55，1对应
% 0.35，naishidu = -10*r2 + 6;




