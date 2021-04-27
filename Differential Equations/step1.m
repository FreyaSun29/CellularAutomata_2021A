function dx=step1(t,x,r1,r2,yanshen,lambda1,lambda2)
    dx=zeros(2,1); 
    dx(1)=x(1)*(r1-lambda1*x(2));
    dx(2)=x(2)*(-r2+lambda2*x(1))+yanshen;
end