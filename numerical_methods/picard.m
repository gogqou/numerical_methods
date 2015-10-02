function res = picard (tridiag, known, guess,x0, range, n, tol)
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
        known = (deltax)^2.*(3*u+10.*u.^3+x.^2)
    end
    if nargin < 1
        tridiag = tridiagnonlinr(n);
    end
    
    dx = range/(n+1);
    temp = ones(n,1);
    xs = x0+dx:dx:x0+range-dx;
    i = 1;
    vector = zeros(1000, n+2);
    vector(1,:) = guess;
    difference = 1;
    while norm(difference)>tol && i< 200
        us=vector(i,:)
        answr = subs(known, {deltax x u}, {dx xs us(2:length(us)-1)})
        answrbound= [0 answr 0]';
        vector(i+1,:) = tridiag\answrbound;
        newguess =vector(i+1,:);
        i = i+1;
        difference = vector(i,:)-vector(i-1,:);
        norm(difference)
       
    end
    i
    x_=[0 xs 1]';
    plot(x_, newguess');
    res=[x_, newguess'];
end


