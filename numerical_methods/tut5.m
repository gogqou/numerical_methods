function res = tut5(ODE,range,n,x0,tol, method)
    clf
    syms deltax x u uprime
    if nargin < 6
        method = @picard;
    end
    
    if nargin <5
        tol = 1e-3;
    end
    if nargin < 4
        x0 = 0;
    end
    if nargin < 3
        n = 100;
    end
    if nargin < 2
        range = 1;
    end
    if nargin < 1
        ODE = 3.*u +10*u.^3 + x.^2 + uprime*0;
    end 
    A=(deltax)^2.*(3*u+10.*u.^3+x.^2);
    points = method(tridiagnonlinr(n), A, zeros(n+2,1),x0,range, n, tol);
    
    plot(points);
    
    
    
    
end