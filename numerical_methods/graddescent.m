function res = graddescent (n, iteratns)
    if nargin < 2
        iteratns =15;
    end
    if nargin < 1
        n = 2;
    end
    guess = [4.5 1.75 .1 9.7]';
    xs = linspace(1,100,20);
    W = eye(length(xs)); %setting weighting all to 1 for now
    p=guess;
    syms x 
    functn_yhat= datapts(n,p);
    dp = .00001;
    func_y_real=datapts(n);
    datapoints=(subs(func_y_real, x, xs))';
    
    for i = 1:iteratns
        fitfunctnpoints=(subs(functn_yhat, x,xs))';
        diff = datapoints- fitfunctnpoints;
        J= jacobian (dp, p, xs,functn_yhat);
%         JtransWJ = J'*W*J;
%         JtransW = J'*W*diff;
        h = J'*W*diff;
        p=p+h;
        functn_yhat = datapts(n,p);
    end
    res=p;
end
function res = jacobian (dp, p, xs,functn_yhat)
syms x
    J = zeros(length(xs), length(p));
    for j = 1:length(p)
        del_p = zeros(length(p),1);
        del_p(j) = dp;
        p_jac = p+del_p;
        functn_yhat_jac = datapts(length(p_jac), p_jac);
        for i = 1:length(xs)
            newpyhat= subs(functn_yhat_jac, x, xs(i));
            oldpyhat =subs(functn_yhat, x, xs(i));
            J(i,j)=(newpyhat - oldpyhat)/dp;
        end
    end
    res = J;
end
function res = datapts (n,p)
    if nargin < 2
       p = [5 2 .2 10]; 
    end
    if nargin < 1
        n= 2;
    end
    syms x
%     func=0;
    func = n/n*p(1) * exp(-x/p(2))+ p(3) *x *exp(-x/p(4));
%         for i = 1:n
%             func=func+sin(p(i)*x);
%         end
    res=func;
end