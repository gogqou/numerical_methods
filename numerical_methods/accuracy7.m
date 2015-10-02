function res = accuracy7(method, K, time, x0, range,u0,N)
syms k x t jay 
    if nargin < 7
        N = 20;
    end
    if nargin < 6
        u0 = [1/2,1];
    end
    if nargin< 5
        range = 1;
    end
    if nargin < 4
        x0 = 0;
    end
    
    if nargin < 3
        time = .05;    
    end
    if nargin< 2
        K = 60;
    end
    if nargin < 1
        method =@galerkin_heat;
    end
    ks = 1: K;
    exact = @exactsolution7;
    dx = range/(K-1);
    xs = x0:dx:x0+range;
    xs = linspace(x0,range);
    for x = 1:length(xs)
        for k = 1: length(ks)
            B(x,:)= exact(ks,xs(x),0);
        end
    end
    u_initial = zeros(length(xs),1);
    u_initial(round(length(xs)* u0(1)),1) = u0(2);
	ak = linsolve(B, u_initial);
     for x = 1:length(xs)
        B(x,:)= exact(ks,xs(x),time);
     end
     exactpts=B*ak;
     exact_pt=exactpts(round(length(xs)* u0(1)));
     for i = 1:N
         N (i) = i;
         
         syms jay x 
         approx=collocation_heat(sin(jay*pi*x), [1/2 1], [0 0], 0,1,i, .0001,.05);
         diff(i)=abs(exact_pt- approx);
     end
    plot(N, diff);
    res =ak;
    
end


function res = exactsolution7(k,x,t)
if nargin < 3
    t = 0;
end
if nargin < 2 
    x = 1;
end
if nargin < 1
    k = 1;
end
res=exp(-k.^2 * pi^2*t).*sin(k*pi*x);
end