function res = accuracytut7 (method, N, dt, x0, range,u0)
syms k x t 
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
        dt = .001;    
    end
    if nargin< 2
        N = 11;
    end
    if nargin < 1
        method =@collocation_ht;
    end
    exactsolutn = exp(-k^2 * pi^2*t)*sin(k*pi*x);
    ks = 1: N;
    exactwithksubbed=subs(exactsolutn, k, ks);
    dx = range/(N+1);
    xs = (x0+dx):dx:(x0+range-dx);
    for i = 1:length(xs)
      exact(i,:)=subs(exactwithksubbed, x, xs(i));
       
    end
    B=subs(exact,t,0);
    u_initial = zeros(1,N);
    u_initial(round(N* u0(1))) = u0(2);
    ak=B/u_initial;
    plot(xs,B*ak);
%     B*ak
    res =ak;


end