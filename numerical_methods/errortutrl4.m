function res=errortutrl4 (range, tolrange, x0, highest)
% clf
    if nargin < 4
        highest= 6;
    end
    if nargin <3
        x0 = 0;
    end
    if nargin < 2
        tolrange = 3;
    end
    if nargin < 1
        range = 1;
    end
    exact = (tut4(range,10^(highest), x0));

        for i = 1:1:tolrange;
            numdiv = 10^i
            factor = 10^highest/numdiv
            x = 1+factor: factor:10^highest-factor
            temp1= tut4(range,numdiv,x0);
            vs = temp1(2:length(temp1)-1);
            length(vs)-length(x);
            real = exact(x)
            dx(i)=range/(numdiv+1);
            diff(i) = norm (real-vs) /sqrt(length(vs));
        end
        clf
    loglog ((dx), (diff), 'o');
     a = log(dx);
     b = log(diff);
%     c = 2*log(dx);
%     hold on
%     plot(a,a, 'r');
%     plot (a,c);
    slope = (b(length(dx)-1)-b(1))/(a(length(dx)-1)-a(1));
    res =[slope (dx); 0 (diff)];

end
    