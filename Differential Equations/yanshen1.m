% 不同温度也影响r2
list_all = zeros(1,101);
% 当其他条件（湿度）不变时，只生长速率（延申）不同的情况下有
for j=0.5:0.05:0.6
    list = [];
    k=1;
    for i=0:0.1:10
        r1=0.00001;
        r2=0.55;
        yanshen=i;
        lambda1=0.00005;
        lambda2=j;
        [t1,x1]=ode45(@(t1,x1)step1(t1,x1,r1,r2,yanshen,lambda1,lambda2),[1 122],[1 0]);
        wood=x1(:,1); fungi=x1(:,2); 
        list(k) = wood(end);
        k = k+1;
    end
    list_all = [list_all;list];
end

plot(0:0.1:10,1-list_all(2,:),"r",0:0.1:10,1-list_all(3,:),"b",0:0.1:10,1-list_all(4,:),"y")
axis([0,10,0,0.6]);

