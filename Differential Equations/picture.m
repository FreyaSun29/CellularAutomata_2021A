clear;clc
list = zeros(21,61);
m = 1;
for j=0.5:0.01:0.7
    k=1;
    for i=0:0.1:6
        r1=0.00001;
        r2=j;
        yanshen=i;
        lambda1=0.00005;
        lambda2=0.5;
        [t1,x1]=ode45(@(t1,x1)step1(t1,x1,r1,r2,yanshen,lambda1,lambda2),[1 122],[1 0]);
        wood2=x1(:,1); fungi2=x1(:,2); 
        list(m,k) = 1-wood2(end);
        k = k+1;
    end
    m=m+1;
end
r2s = 0.5:0.01:0.7;
x = -10*r2s + 6;
y = 0:0.1:6;
z = list;

subplot(1,2,1);
mesh(y,x,z)
xlabel('Growth Rate');  ylabel('Moisture Tolerance');  zlabel('Decomposition Rate');
title('Breakdown Model(Mesh)');

subplot(1,2,2);
surf(y,x,z)
shading interp;
xlabel('Growth Rate');  ylabel('Moisture Tolerance');  zlabel('Decomposition Rate');
title('Breakdown Model(surf)');





