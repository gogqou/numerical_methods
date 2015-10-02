function res = tridiagnonlinr(n)
    if nargin < 1
        n =6;
    end
      
%      xs = x0+dx: dx: x0+range-dx;
    temp = ones(n,1)';
    main =[-2*temp];
    low = ones(n-1,1)';
    up = [1 temp(3:end)];
   
    res = diag([0 low 0],-1)+diag([1 main 1]) + diag([0 up 0], 1);
%     res = [up 0; main; 0 low; known];
%     res = [up 0; main; 0 low];
end