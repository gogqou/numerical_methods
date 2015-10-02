function res = finitediff(N,dt,u0, x0, range, bound,time, c)
syms 
    if nargin < 8
        c = 1;
    end
    if nargin < 7
        time = .01; % how long to run: determines highest value of m
    end

    if nargin < 6
        bound = [0 0]; %in heat eqtn case--the ends always have temp = 0
    end    
    if nargin < 5
        range =1; % normalized length of rod: 0 to 1
    end
    if nargin < 4
        x0 = 0;
    end
    if nargin <3
        u0=[1/2 1]; %   u0(1) is the location along the localized length [0 1] and u0(2) is the actual value
    end
    if nargin < 2
        dt = .0001;
    end    
    if nargin<1
        N=40;
    end
    %% u_m = (I + dt*A)^m * u0
    dx = range/(N+1);
    u_init = zeros(1,N); %sets up initial values generally--all zero
    u_init(round(u0(1) * N)) = u0(2); % sets the initial value with givens; u0(1) * N gives index
    A_constantpart = c/(dx^2);
    maindiag= (-2* ones(N,1));
    sidediag= ones(N-1,1);
    A_matrix = diag(maindiag)+diag(sidediag,-1)+diag(sidediag,1);
    A= A_constantpart* A_matrix;
    IplusdttimesA = eye(N)+ dt* A;
    for m = 1:(time/dt)
       time (m) = dt*m;
       u(m,:)=[bound(1) u_init*(IplusdttimesA)^m bound(2)];
    end
    xs = x0: dx: (x0+range);
    surf(xs, time, u);
    set(0, 'defaultaxesfontsize', 18)
    xlabel('X', 'FontSize', 20);
    ylabel('Time(s)', 'FontSize', 20);
    zlabel ('Temperature', 'FontSize', 20);
    res = IplusdttimesA;


end   
