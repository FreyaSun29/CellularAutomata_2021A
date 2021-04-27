function dx=step1(t,x,rate,humidity_resistance)
    r1=0.00001;
    naishidu = humidity_resistance;
    r2=-(naishidu-6)/10;
    yanshen=rate;
    lambda1=0.00005;
    lambda2=0.5;
    dx=zeros(2,1); 
    dx(1)=x(1)*(r1-lambda1*x(2));
    dx(2)=x(2)*(-r2+lambda2*x(1))+yanshen;
end