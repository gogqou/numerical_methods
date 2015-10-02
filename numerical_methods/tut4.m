function res= tut4(range, n,x0)
if nargin < 3
    x0 = 0;
end
if nargin < 2
    n = 100;
end
if nargin < 3
    range = 1;
end
    dx = range/(n-1);
    matrices =diagmatrixbound ([1 0 1], [1 0 0], n);
    bottom = matrices (3,:);
    main = matrices(2,:);
    up = matrices (1,:);
    known = matrices (4,:);
    vs=tridiagonal (bottom(2:end), main, up(1:end-1), known);
    xs= [x0:range/(n-1):x0+range];
    plot (xs,vs);
    res =[vs];
end
