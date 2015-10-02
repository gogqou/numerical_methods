function res = collocation (ODE, phi, bound, x0, range, n)
    syms x u l uprime udprime phi phiprime phidprime

    if nargin < 6
        n = 20;
    end
    if nargin <5
        range = 1;
    end
    if nargin < 4
        x0 = 0;
    end
    if nargin < 3
        bound = [0 0];
    end
    if nargin < 2
%         phi = (x-1)*(1+x)^l;
        phi = sin(l*pi*x);
%          phi = cos(l*pi*x);
    end
    if nargin < 1
        %%udprime%%
        ODE = 4*exp(4*(2*x-1));
%         ODE = exp(4*x);
    end
    
   %%general setup: full matrix of phi(x) * As = the known; 
   %%evaluate the basis function matrix, evaluate the "known" matrix and
   %%solve for As--coefficients?
    ls = linspace(1,n,n)';
    xs = linspace(x0+range/n,x0+range-range/n,n);
%     xs = linspace(x0, x0+range, n+2);
%     xs = xs(2:n+1);
    phiprime = diff(phi,x);
    phidprime = diff (phiprime,x); %phidoubleprime
    a=subs(phidprime, x, xs); % step wise substitution, first substitute xs
    phimatrix=subs(a,l,ls)'; % then substitute ls
    known = subs(ODE, x, xs'); % calculate known vector
    As = linsolve(phimatrix,known); % calculate coefficients
%     As =  phimatrix\known
    solutn= As*phi;
        for i = 1:n;
            b(i)=subs(solutn(i), l, i);
        end
    
%     Js = 1:n;
%     fcn = matlabFunction(phi);
%     
%     for i =1:n
%         Y(i) = fcn(Js, xs(i))*As;
%     end
    
    c = sum(b);
    eqtn= matlabFunction(c);
    pts= eqtn(xs);
    xs_ = [x0 xs x0+range];
    pts_ = [bound(1) pts bound(2)];
    plot(xs_, pts_);
    res=eqtn;
       
end
        