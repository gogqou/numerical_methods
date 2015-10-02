function res = galerkin (ODE, basis, bound, x0, range, n)
syms x u l uprime udprime basis phi phiprime 

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
%         basis = x*(1-x)^l;
        basis = 2^(1/2)/(l*pi)*sin(l*pi*x);
%         basis = sin(l*pi*x);
%          basis = cos(l*pi*x);
    end
    if nargin < 1
        %%udprime%%
        ODE = 4*exp(4*(2*x-1));
%         ODE = exp(4*x);
    end
%% setup: sum j = 1 to N times coefficients times integral from 0 to 1 
%% -integral of (phi_i (x) * phi_j(x)) for i = 1, ..., N
%% sets up a diagonal matrix 
%%multiply by the coefficient vector a_j
%%multiply the matrix by the vector a_j, this equals the integral from 0 to
%%1 of the second derivative of the function times the basis function
%%integrand * a_j = d(x) * phi_i(x)
    phi= basis;
    phiprime = diff(phi, x);
    integrand = phiprime*phiprime; % phiprime * phiprime
    integrand2 = ODE*phi;
    diags = zeros (1,n);
    answr = zeros(1,n);
    for i = 1:n
        b = subs(integrand, l, i);
        c = matlabFunction(b);
        d = subs(integrand2, l, i);
        f = matlabFunction(d);
        diags(i) = quad(c,0,1);
        answr(i) = quad(f,0,1);
    end
    phi2 = -diag(diags) % makes the diagonal matrix--only true in this setup
    a=linsolve(phi2,answr'); %divide phi2 into answer vector
    basis_matrix_standin=a*basis;
    
    for i = 1:n
        basismatrix(i)= subs(basis_matrix_standin(i), l, i);
    end
    expressionforU=sum(basismatrix);
    xs = x0+range/n :range/n: x0+range-range/n;
    eqtn = matlabFunction(expressionforU);
    pts = eqtn(xs);
    xs_= [x0 xs x0+range];
    pts_= [bound(1) pts bound(2)];
    plot(xs_,pts_);
    res=expressionforU;
    
    
     
    
    
    
    
    
    
end