function res = nonlinnewton(tridiag, known, guess, x0, range, n, tol)
    syms deltax x u uprime
    if nargin < 7
        tol = 1e-3;
    end
    if nargin < 6
        n = 100;
    end
    if nargin <5
        range = 1;
    end
    if nargin < 4
        x0 = 0;
    end
    if nargin < 3
        guess = zeros(n+2,1);
    end
    if nargin < 2
        known = (deltax)^2.*(3*u+10.*u.^3+x.^2);
    end
    if nargin < 1
        tridiag = tridiagnonlinr(n);
    end
    dx = range/(n+1);
    temp = ones(n,1);
    dxs= dx.*temp;
    xs = x0+dx:dx:x0+range-dx;
    vars=symvar(known);
    Dknown = diff(known,vars)
    i = 1;
    us =guess;
    delta = 1;
    while norm(delta)>tol
        %%F is Au - H(u) in the setup Au = H(u); DF is derivative of F%% 
        U=us(2:length(us)-1);
        F=subs(known, {deltax x u}, {dx xs U'});
%         F = tridiag*us -subs(known, {deltax x u}, {dxs, xs,us});
        random = [U U U];
        %% test stuff to figure out why the matrices aren't the right size
        r=ones(length(U),1)
        subs([u^3 30*u u+x^2], {x u}, {r U})
        subs(Dknown, {u}, {U})
        %%
        DF = tridiag - subs(Dknown, {deltax u}, {dxs U})
        delta= linsolve(DF, -F);
        us = us + delta';
    end
    
    res = us;
end
    
    