function res = bessel_()
n = 10000;
s = 10*10^-3;
w=15*10^(-6)*10^3;
e0=8.85*10^-14;
er=2.8;
eair=1;
L= 600*10^-3;
N= 25;
for i= 1:n
    term1 = 1/(2*i-1);
    top=(2*i-1)*pi*s;
    bottom = 2*(s+w);
    a= besselj(0, top/bottom);
    b=a^2;
    c(i)=term1*b;
end
res = L*N*4*e0*eair/pi*sum(c);
    
