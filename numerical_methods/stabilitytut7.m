function res = stabilitytut7 (method, N, dt)
    if nargin < 3
        dt = [1e-4 .01];
    end
    if nargin< 2
        N =[2 50];
    end
    if nargin < 1
        method = @finitediff;
    end

    dts = dt(1):1e-3: dt(2);
    ns = N(1) : N(2);
    EVmax = zeros(length(ns), length(dts));
    for nindex = 1:length(ns)
        for tindex = 1:length(dts) 
            tstepsize=dts(tindex);
            EV= eig(method(ns(nindex), tstepsize));
            EVmax(nindex, tindex)=max(EV);
        end
    end
%     size(dts)
%     size(ns)
%     size(EVmax)
    surf(dts,ns,EVmax);
    set(0, 'defaultaxesfontsize', 18)
    xlabel('Stepsize', 'FontSize', 20);
    ylabel('N', 'FontSize', 20);
    zlabel ('Eigenvalue', 'FontSize', 20);
    
    
end


